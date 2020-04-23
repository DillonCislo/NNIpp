/*
 * Copyright (C) 2020 Dillon Cislo
 *
 * This file is part of NNI++.
 *
 * NNI++ is free software: you can redistribute it and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will by useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTIBILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program.
 * If not, see <http://www.gnu.org/licenses/>
 *
 */

#include "NaturalNeighborInterpolant.h"

#include <cassert>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

#include <igl/edges.h>
#include <igl/adjacency_list.h>
#include <igl/boundary_loop.h>

#include "ghostPointsCircle.h"
#include "ghostPointsEdge.h"
#include "../General/delaunayTriangulation.h"
#include "../General/convexHull.h"
#include "../General/inPolygon.h"
#include "../General/triCircumcenter.h"
#include "../General/sortPolygonCCW.h"
#include "../General/polyArea.h"
#include "delta.h"
#include "gradientsDirect.h"
#include "hessianDirect.h"

// ======================================================================================
// CONSTRUCTOR FUNCTIONS
// ======================================================================================

///
/// Default constructor
///
template <typename Scalar>
NNIpp::NaturalNeighborInterpolant<Scalar>::NaturalNeighborInterpolant(
    const Vector &Xp, const Vector &Yp,
    const Matrix &Vp, const NNIParam<Scalar> &param ) {

  // Checck for valid parameters
  param.checkParam();

  // Set display parameters
  this->m_dispGrad = param.dispGrad;
  this->m_dispInterp = param.dispInterp;

  // Generate the ghost points for extrapolation
  Vector GPx;
  Vector GPy;
  this->generateGhostPoints( Xp, Yp, param, GPx, GPy );

  // Combine data points and ghost points into the extended points list
  int numAllPoints = param.numPoints + GPx.size();
  this->m_Points =
    Eigen::Matrix<Scalar, Eigen::Dynamic, 2>::Zero( numAllPoints, 2 );
  this->m_Points.col(0) << Xp, GPx;
  this->m_Points.col(1) << Yp, GPy;
  
  // Construct the extended triangulation
  Vector Px = this->m_Points.col(0);
  Vector Py = this->m_Points.col(1);
  this->m_Faces = Eigen::Matrix<int, Eigen::Dynamic, 3>::Zero(1,3);
  NNIpp::delaunayTriangulation( Px, Py, this->m_Faces );

  // Construct the extended triangulation edge list
  this->m_Edges = Eigen::Matrix<int, Eigen::Dynamic, 2>::Zero(1,2);
  igl::edges( this->m_Faces, this->m_Edges );

  // Calculate the circumcenter of each face of the extended triangulation
  int numFaces = this->m_Faces.rows();
  Eigen::Matrix<Scalar, Eigen::Dynamic, 3> FVx( numFaces, 3 );
  Eigen::Matrix<Scalar, Eigen::Dynamic, 3> FVy( numFaces, 3 );

  for( int i = 0; i < numFaces; i++ ) {

    int vi = this->m_Faces(i,0);
    int vj = this->m_Faces(i,1);
    int vk = this->m_Faces(i,2);

    FVx.row(i) << Px(vi), Px(vj), Px(vk);
    FVy.row(i) << Py(vi), Py(vj), Py(vk);

  }

  NNIpp::triCircumcenter( FVx, FVy, this->m_FCC );
  this->m_FaceVertexX = FVx;
  this->m_FaceVertexY = FVy;

  // Calculate the convex hull of the extended triangulation
  this->m_ConvexHull = Eigen::Matrix<Scalar, Eigen::Dynamic, 2>::Zero(1,2);
  NNIpp::convexHull( Px, Py, this->m_ConvexHull );

  // Calculate the 'delta' parameters for each face of the extended triangulation
  NNIpp::delta( FVx, FVy, this->m_fDelta );

  // Precompute data for calculating the 'gamma' parameter of each face
  this->precomputeGamma();

  // Update the basic values matrix
  this->m_Values = Matrix::Zero( numAllPoints, param.numVals );
  for( int i = 0; i < param.numPoints; i++ ) {
    for( int j = 0; j < param.numVals; j++ ) {
      this->m_Values(i,j) = Vp(i,j);
    }
  }

  // Generate the gradients of the input points
  this->generateGradients( param );

  // Set the values/gradients of the ghost points
  this->ghostPointValueHandling( param );

};

///
/// Generate ghost points
///
template <typename Scalar>
NNI_INLINE void NNIpp::NaturalNeighborInterpolant<Scalar>::generateGhostPoints(
    const Vector &Xp, const Vector &Yp,
    const NNIParam<Scalar> &param,
    Vector &GPx, Vector &GPy ) {

  // Choose the appropriate ghost point generation algorithm
  switch( param.ghostMethod ) {

    case NNI_GHOST_POINTS_CUSTOM :

      GPx = param.customGhostPoints.col(0);
      GPy = param.customGhostPoints.col(1);
      break;

    case NNI_GHOST_POINTS_CIRCLE :

      NNIpp::ghostPointsCircle( Xp, Yp, param, GPx, GPy );
      break;

    case NNI_GHOST_POINTS_EDGE :

      NNIpp::ghostPointsEdge( Xp, Yp, param, GPx, GPy );
      break;

  }

  // Check that all of the data points lie within the convex hull of the ghost points
  int GPn = GPx.size();
  Eigen::Matrix<Scalar, Eigen::Dynamic, 2> GPCH( GPn, 2 );
  NNIpp::convexHull( GPx, GPy, GPCH );

  Eigen::Array<bool, Eigen::Dynamic, 1> inPoly(GPn, 1);
  NNIpp::inPolygon( GPCH, Xp, Yp, inPoly );

  for( int i = 0; i < GPn; i++ ) {
    if (!inPoly(i)) {
      std::runtime_error("Invalid ghost point convex hull");
    }
  }

};

///
/// Generate discrete gradients
///
template <typename Scalar>
NNI_INLINE void NNIpp::NaturalNeighborInterpolant<Scalar>::generateGradients(
    const NNIParam<Scalar> &param ) {

  // Allocate memory for derivative arrays
  int numAllPoints = this->m_Values.rows();
  this->m_DataGradX = Matrix::Zero( numAllPoints, param.numVals );
  this->m_DataGradY = Matrix::Zero( numAllPoints, param.numVals );
  this->m_DataHessXX = Matrix::Zero( numAllPoints, param.numVals );
  this->m_DataHessXY = Matrix::Zero( numAllPoints, param.numVals );
  this->m_DataHessYY = Matrix::Zero( numAllPoints, param.numVals );

  if ( param.gradientSupplied ) {

    // Set user supplied gradient data
    for( int i = 0; i < param.numPoints; i++ ) {
      for( int j = 0; j < param.numVals; j++ ) {

        this->m_DataGradX(i,j) = param.DataGradX(i,j);
        this->m_DataGradY(i,j) = param.DataGradY(i,j);

      }
    }

    if ( param.hessianSupplied ) {

      // Set user supplied Hessian data
      for( int i = 0; i < param.numPoints; i++ ) {
        for( int j = 0; j < param.numVals; j++ ) {

          this->m_DataHessXX(i,j) = param.DataHessXX(i,j);
          this->m_DataHessXY(i,j) = param.DataHessXY(i,j);
          this->m_DataHessYY(i,j) = param.DataHessYY(i,j);

        }
      }

    } else {

      // Generate discrete Hessian data
      switch( param.gradType ) {

        case NNI_GRADIENTS_DIRECT :

          NNIpp::hessianDirect( this->m_Points, this->m_Faces,
              this->m_Values, param, this->m_DataGradX, this->m_DataGradY,
              this->m_DataHessXX, this->m_DataHessXY, this->m_DataHessYY );
          
          break;

        case NNI_GRADIENTS_ITER :

          std::runtime_error("This functionality is not available yet");
          break;

      }

    }

  } else {

    // Generate discrete gradient and Hessian data
    switch( param.gradType ) {

      case NNI_GRADIENTS_DIRECT :

        NNIpp::gradientsDirect( this->m_Points, this->m_Faces,
            this->m_Values, param,
            this->m_DataGradX, this->m_DataGradY,
            this->m_DataHessXX, this->m_DataHessXY, this->m_DataHessYY );

        break;

      case NNI_GRADIENTS_ITER :

        std::runtime_error("This functionality is not available yet");
        break;

    }

  }

};

///
/// Generate values/derivatives for ghost points
///
template <typename Scalar>
NNI_INLINE void NNIpp::NaturalNeighborInterpolant<Scalar>::ghostPointValueHandling(
    const NNIParam<Scalar> &param ) {

  // The number of points in the extended triangulation
  int numAllPoints = this->m_Points.rows();

  // The number of input data points
  int numPoints = param.numPoints;

  // The number of input value components
  int numVals = param.numVals;

  // The number of ghost points
  int numGP = numAllPoints - numPoints;

  // Construct the vertex adjacency list of the extended triangulation
  std::vector<std::vector<int> > VA( numAllPoints );
  igl::adjacency_list( this->m_Faces, VA );

  // Iterate over all ghost points
  for( int i = numPoints; i < numAllPoints; i++ ) {

    // The vertices adjacent to the current ghost point
    std::vector<int> curVA = VA[i];

    // Remove other ghost points from consideration
    std::vector<int>::iterator it;
    it = std::remove_if( curVA.begin(), curVA.end(),
        std::bind2nd( std::greater<int>(), (numPoints-1) ) );
    curVA.erase( it, curVA.end() );

    // The number of non-ghost points attached to the current ghost point
    int numVA = curVA.size();

    // Calculate edge vectors and inverse edge lenths
    Vector X( numVA );
    Vector Y( numVA );
    Vector W( numVA );

    for( int k = 0; k < numVA; k++ ) {

      X(k) = this->m_Points(curVA[k], 0) - this->m_Points(i, 0);
      Y(k) = this->m_Points(curVA[k], 1) - this->m_Points(i, 1);
      W(k) = 1 / std::sqrt( X(k)*X(k) + Y(k)*Y(k) );

    }

    // Edge weights are normalized inverse edge lengths
    W = ( W.array() / W.sum() ).matrix();

    // Iterate over all input value components
    for( int j = 0; j < numVals; j++ ) {

      Scalar Vij = Scalar(0.0);
      Scalar DXij = Scalar(0.0);
      Scalar DYij = Scalar(0.0);
      Scalar DXXij = Scalar(0.0);
      Scalar DXYij = Scalar(0.0);
      Scalar DYYij = Scalar(0.0);

      // Iterate over the vertices adjacent to the current ghost point
      for( int k = 0; k < numVA; k++ ) {

        // Calculate contribution to average value
        Vij += W(k) * this->m_Values(curVA[k], j);

        // Calculate contribution to average gradients
        DXij += W(k) * ( this->m_DataGradX(curVA[k], j) +
            this->m_DataHessXX(curVA[k], j) * X(k) +
            this->m_DataHessXY(curVA[k], j) * Y(k) );

        DYij += W(k) * ( this->m_DataGradY(curVA[k], j) +
            this->m_DataHessXY(curVA[k], j) * X(k) +
            this->m_DataHessYY(curVA[k], j) * Y(k) );

        // Calculate contribution to average Hessians
        DXXij += W(k) * this->m_DataHessXX(curVA[k], j);
        DXYij += W(k) * this->m_DataHessXY(curVA[k], j);
        DYYij += W(k) * this->m_DataHessYY(curVA[k], j);
        
      }

      // Update ghost points values/gradients
      this->m_Values(i,j) = Vij;
      this->m_DataGradX(i,j) = DXij;
      this->m_DataGradY(i,j) = DYij;
      this->m_DataHessXX(i,j) = DXXij;
      this->m_DataHessXY(i,j) = DXYij;
      this->m_DataHessYY(i,j) = DYYij;

    }

  }

};

///
/// Precompute data needed to calculate the 'Gamma' parameter
/// on each face of the extended triangulation
///
template <typename Scalar>
NNI_INLINE void NNIpp::NaturalNeighborInterpolant<Scalar>::precomputeGamma() {

  ArrayVec x1 = this->m_FaceVertexX.col(0);
  ArrayVec x2 = this->m_FaceVertexX.col(1);
  ArrayVec x3 = this->m_FaceVertexX.col(2);

  ArrayVec y1 = this->m_FaceVertexY.col(0);
  ArrayVec y2 = this->m_FaceVertexY.col(1);
  ArrayVec y3 = this->m_FaceVertexY.col(2);

  ArrayVec r1S = x1*x1 + y1*y1;
  ArrayVec r2S = x2*x2 + y2*y2;
  ArrayVec r3S = x3*x3 + y3*y3;

  ArrayVec G1 = r1S * ( x2 * y3 - y2 * x3 );
  G1 += r2S * ( x3 * y1 - y3 * x1 );
  G1 += r3S * ( x1 * y2 - y1 * x2 );

  ArrayVec G2 = r1S * (y2-y3) + r2S * (y3-y1) + r3S * (y1-y2);
  ArrayVec G3 = r1S * (x3-x2) + r2S * (x1-x3) + r3S * (x2-x1);
  ArrayVec G4 = x1 * (y3-y2) + x2 * (y1-y3) + x3 * (y2-y1);

  Eigen::Array<Scalar, Eigen::Dynamic, 4> G( this->m_Faces.rows(), 4 );
  G << G1, G2, G3, G4;

  this->m_fGamma = G;

};

///
/// Calculate the 'Gamma' parameter from (Hiyoshi, 2008). Used to
/// calculate the natural neighbor coordinates of a query point
/// Gamma(v1, v2, v3, v4) does not have a simply geometric interpretation
///
template <typename Scalar>
NNI_INLINE void NNIpp::NaturalNeighborInterpolant<Scalar>::gamma(
    const Scalar X, const Scalar Y, ArrayVec &G ) {

  int numF = this->m_Faces.rows();

  ArrayVec XX = ArrayVec::Constant( numF, 1, X );
  ArrayVec YY = ArrayVec::Constant( numF, 1, Y );

  G = this->m_fGamma.col(0) +
    this->m_fGamma.col(1) * XX +
    this->m_fGamma.col(2) * YY +
    this->m_fGamma.col(3) * ( XX*XX + YY*YY );

};

// ======================================================================================
// NATURAL NEIGHBOR COORDINATE FUNTIONS
// ======================================================================================

///
/// Calculate the natural neighbor coordinates of a set of query points
/// Also extract the sufficient geometric data to calculate the derivative
/// of the natural neighbor coordinates with respect to the query point
/// coordinates
///
/// NOTE: IT IS ASSUMED THAT THE QUERY POINTS LIE WITHIN THE CONVEX HULL
/// OF THE EXTENDED TRIANGULATION
///
template <typename Scalar>
void NNIpp::NaturalNeighborInterpolant<Scalar>::naturalNeighborCoordinates(
    const Vector &Xq, const Vector &Yq,
    std::vector<Vector> &u, std::vector<Eigen::VectorXi> &uIDx,
    std::vector<Eigen::Matrix<Scalar, Eigen::Dynamic, 2> > &uVC,
    Vector &uA ) {

  typedef Eigen::Array<bool, Eigen::Dynamic, 1> BoolVec;

  int numQ = Xq.size(); // The number of query points
  int numP = this->m_Points.rows(); // The number of vertices
  int numF = this->m_Faces.rows(); // The number of faces

  // Iterate over query points
  for( int i = 0; i < numQ; i++ ) {

    //-----------------------------------------------------------------------------------
    // Handle the case that a query point is equal to a control point
    //-----------------------------------------------------------------------------------
    
    // Find the shortest distance between the current query point and ANY control point
    Scalar dx = this->m_Points(0,0) - Xq(i);
    Scalar dy = this->m_Points(0,1) - Yq(i);
    Scalar minDist = std::sqrt( dx*dx + dy*dy );
    int minID = 0;

    for( int j = 1; j < numP; j++ ) {

      dx = this->m_Points(j,0) - Xq(i);
      dy = this->m_Points(j,1) - Yq(i);
      Scalar curDist = std::sqrt( dx*dx + dy*dy );

      if ( curDist < minDist ) {
        minDist = curDist;
        minID = j;
      }

    }

    if ( minDist < Scalar(1.0e-12) ) {

      u[i] = Vector::Constant(1, 1, Scalar(1.0));
      uIDx[i] = Eigen::VectorXi::Constant(1, 1, minID);
      uVC[i] = Eigen::Matrix<Scalar, Eigen::Dynamic, 2>::Constant(1, 2, Scalar(0.0));
      uA(i) = Scalar(0.0);

    }

    //-----------------------------------------------------------------------------------
    // Extract the natural neighbors of the current query point
    //-----------------------------------------------------------------------------------
    
    // The 'Gamma' parameter of each face with respect to the query point
    ArrayVec fGamma( numF, 1 );
    this->gamma( Xq(i), Yq(i), fGamma );

    // The 'Xi' parameter of each face with respect to the query point
    // The faces with Xi > 0 contain the query point in their circumcircle
    ArrayVec fXi = fGamma / this->m_fDelta;

    // The union of these faces is a star-shaped 'natural neighbor polygon'
    // The vertices of this polygon are the natural neighbors of the query point
    
    // Extract the number of faces in the natural neighbor polygon
    int numInFace = 0;
    for( int j = 0; j < numF; j++ ) {
      if ( fXi(j) > Scalar(0.0) ) {
        numInFace++;
      }
    }
    
    // The connectivity of those faces
    Eigen::Matrix<int, Eigen::Dynamic, 3> nnPoly_F( numInFace, 3 );

    // The circumcenters of those faces
    Eigen::Matrix<Scalar, Eigen::Dynamic, 2> P_CC( numInFace, 2 );

    int fCount = 0;
    for( int j = 0; j < numF; j++ ) {
      if ( fXi(j) > Scalar(0.0) ) {
        nnPoly_F.row(fCount) = this->m_Faces.row(j);
        P_CC.row(fCount) = this->m_FCC.row(j);
        fCount++;
      }
    }

    // The (CCW sorted) natural neighbors of the current query point
    Eigen::VectorXi quIDx(1);
    igl::boundary_loop( nnPoly_F, quIDx );

    // The number of natural neighbors
    int numNN = quIDx.size();

    //-----------------------------------------------------------------------------------
    // Extract the virtual Voronoi cell of the current query point
    //-----------------------------------------------------------------------------------
    
    // Extract the (x,y)-coordinates of the triangles of the virtual Delaunay
    // triangulation containing the query point
    Eigen::Matrix<Scalar, Eigen::Dynamic, 3> nnTri_X( numNN, 3 );
    Eigen::Matrix<Scalar, Eigen::Dynamic, 3> nnTri_Y( numNN, 3 );

    for( int j = 0; j < numNN; j++ ) {

      int k = (j+1) % numNN;
      nnTri_X.row(j) << Xq(i), this->m_Points(quIDx(j), 0), this->m_Points(quIDx(k), 0);
      nnTri_Y.row(j) << Yq(i), this->m_Points(quIDx(j), 1), this->m_Points(quIDx(k), 1);

    }

    // The circumcenters of those triangles are the vertices of the virtual Voronoi
    // cell of the query point
    Eigen::Matrix<Scalar, Eigen::Dynamic, 2> quVC( numNN, 2 );
    NNIpp::triCircumcenter( nnTri_X, nnTri_Y, quVC );

    //-----------------------------------------------------------------------------------
    // Calculate Sibson's natural neighbor coordinates for each natural neighbor
    //-----------------------------------------------------------------------------------
    
    // The total area of the virtual Voronoi cell
    Scalar qA = Scalar(0.0);

    // Vector of Sibson's coordinates
    Vector qu = Vector::Zero( numNN, 1 );

    int jPrev = numNN - 1;
    for( int j = 0; j < numNN; j++ ) {

      // Find the number of faces in the natural neighbor polygon
      // that contain the current natural neighbor as a vertex
      int numNNinF = 0;
      for( int k = 0; k < numInFace; k++ ) {

        bool addFace = (nnPoly_F(k,0) == quIDx(j)) ||
          (nnPoly_F(k,1) == quIDx(j)) ||
          (nnPoly_F(k,2) == quIDx(j));

        if ( addFace ) {
          numNNinF++;
        }

      }
      
      // Extract the vertices of the intersection of the virtual Voronoi cell
      // of the query point and the Voronoi cell of the current natural neighbor
      Eigen::Matrix<Scalar, Eigen::Dynamic, 2> Pj_V( numNNinF+2, 2 );
      int vvCount = 0;
      for( int k = 0; k < numInFace; k++ ) {

        bool addFace = (nnPoly_F(k,0) == quIDx(j)) ||
          (nnPoly_F(k,1) == quIDx(j)) ||
          (nnPoly_F(k,2) == quIDx(j));

        if ( addFace ) {
          Pj_V.row(vvCount) = P_CC.row(k);
          vvCount++;
        }

      }

      Pj_V.row(numNNinF) = quVC.row(j);
      Pj_V.row(numNNinF+1) = quVC.row( jPrev );

      jPrev = j;

      // CCW sort the vertices of this polygon
      Eigen::Matrix<Scalar, Eigen::Dynamic, 2> Pj_V_CCW( numNNinF+2, 2 );
      Eigen::VectorXi order( numNNinF+2 );
      NNIpp::sortPolygonCCW( Pj_V, Pj_V_CCW, order );

      // Find the area of this polygon
      Scalar Pj_A = NNIpp::polyArea( Pj_V_CCW );

      qA += Pj_A;
      qu(j) = Pj_A;

    }

    // Normalize the coordinates
    qu = ( qu.array() / qA ).matrix();

    // Update the output parameters
    u[i] = qu;
    uIDx[i] = quIDx;
    uVC[i] = quVC;
    uA(i) = qA;

  }

};
      

//TODO: Add explicit template instantiation
#ifdef NNI_STATIC_LIBRARY
#endif
