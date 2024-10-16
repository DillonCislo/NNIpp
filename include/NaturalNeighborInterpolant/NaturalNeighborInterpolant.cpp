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
#include <list>
#include <omp.h>

#include <igl/edges.h>
#include <igl/adjacency_list.h>
#include <igl/boundary_loop.h>
#include <igl/oriented_facets.h>
#include <igl/triangle_triangle_adjacency.h>

#include "ghostPointsCircle.h"
#include "ghostPointsEdge.h"
#include "../General/delaunayTriangulation.h"
#include "../General/convexHull.h"
#include "../General/inPolygon.h"
#include "../General/triCircumcenter.h"
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

  // Check for valid parameters
  param.checkParam();

  // Set display parameters
  this->m_dispGrad = param.dispGrad;
  this->m_dispInterp = param.dispInterp;

  // Generate the ghost points for extrapolation
  Vector GPx = Vector::Zero(param.GPn, 1);
  Vector GPy = Vector::Zero(param.GPn, 1);
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

  // Construct the oriented edge list of the extended triangulation
  this->m_OrientedEdges = Eigen::Matrix<int, Eigen::Dynamic, 2>::Zero(1,2);
  igl::oriented_facets(this->m_Faces, this->m_OrientedEdges);

  // Construct face adjacency list
  this->m_FaceNeighbors = Eigen::Matrix<int, Eigen::Dynamic, 3>::Zero(1,3);
  Eigen::Matrix<int, Eigen::Dynamic, 3> TTI(1,3);
  igl::triangle_triangle_adjacency( this->m_Faces, this->m_FaceNeighbors, TTI );

  Eigen::PermutationMatrix<3,3> perm(3);
  perm.indices() = Eigen::Vector3i(1,2,0);
  this->m_FaceNeighbors = (this->m_FaceNeighbors * perm).eval();

  // Construct the AABB tree for the extended triangulation
  Matrix tmpP = this->m_Points;
  Eigen::MatrixXi tmpF = this->m_Faces;
  this->m_Tree.init( tmpP, tmpF );

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

  if ( !inPoly.array().all() )
    std::runtime_error("Invalid ghost point convex hull");

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
    curVA.erase(
      std::remove_if(curVA.begin(), curVA.end(),
        [numPoints](int idx) { return idx >= numPoints; }),
      curVA.end());

    /* OLD CODE (deprecated in C++11 and removed in C++17)
    std::vector<int>::iterator it;
    it = std::remove_if( curVA.begin(), curVA.end(),
        std::bind2nd( std::greater<int>(), (numPoints-1) ) );
    curVA.erase( it, curVA.end() );
    */

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
    const Scalar X, const Scalar Y, ArrayVec &G ) const {

  int numF = this->m_Faces.rows();

  ArrayVec XX = ArrayVec::Constant( numF, 1, X );
  ArrayVec YY = ArrayVec::Constant( numF, 1, Y );

  G = this->m_fGamma.col(0) +
    this->m_fGamma.col(1) * XX +
    this->m_fGamma.col(2) * YY +
    this->m_fGamma.col(3) * ( XX*XX + YY*YY );

};

///
/// Determine if a query point lies within the circumcircle
/// of a face of the extended triangulation
///
template <typename Scalar>
NNI_INLINE bool NNIpp::NaturalNeighborInterpolant<Scalar>::inCircle(
    const Scalar Xq, const Scalar Yq, const int FID ) const {

  Scalar gamma = this->m_fGamma(FID, 0) +
    this->m_fGamma(FID, 1) * Xq +
    this->m_fGamma(FID, 2) * Yq +
    this->m_fGamma(FID, 3) * ( Xq*Xq + Yq*Yq );

  Scalar xi = gamma / this->m_fDelta(FID);

  // A query point will lie within the circumcircle of
  // a face if the 'xi' parameter is greater than zero
  return ( xi > Scalar(0.0) );

};

///
/// Determine whether each point in a list of query points lies
/// within the faces of the extended triangulation
///
template <typename Scalar>
NNI_INLINE void NNIpp::NaturalNeighborInterpolant<Scalar>::inElement(
    const Vector &Xq, const Vector &Yq, Eigen::VectorXi &I ) const {

  // The number of query points
  const int numQ = Xq.size();

  I.setConstant( numQ, 1, -1 );

  // Temporary copies are inefficient but required due to
  // poor templating in the 'AABB' class implemented in libigl
  Matrix tmpP = this->m_Points;
  Eigen::MatrixXi tmpF = this->m_Faces;

  #pragma omp parallel for if (numQ > 10000)
  for( int i = 0; i < numQ; i++ ) {

    Eigen::Matrix<Scalar, 1, Eigen::Dynamic> Qrow(1,2);
    Qrow(0) = Xq(i);
    Qrow(1) = Yq(i);

    const auto R = this->m_Tree.find( tmpP, tmpF, Qrow, true );

    if (!R.empty()) {
      I(i) = R[0];
    }

  }

  // Check for out-of-bounds query points
  for( int i = 0; i < numQ; i++ ) {
    if (I(i) == -1) {
      throw std::runtime_error("Query points lie outside convex hull of control points!");
    }
  }

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
template <typename Scalar>
void NNIpp::NaturalNeighborInterpolant<Scalar>::naturalNeighborCoordinates(
    const Vector &Xq, const Vector &Yq,
    std::vector<Vector> &u, std::vector<Eigen::VectorXi> &uIDx,
    std::vector<Eigen::Matrix<Scalar, Eigen::Dynamic, 2> > &uVC,
    Vector &uA ) const {

  typedef Eigen::Array<bool, Eigen::Dynamic, 1> VectorXb;
  typedef Eigen::Matrix<Scalar, 1, 2> RowVec2d;

  int numQ = Xq.size(); // The number of query points
  int numP = this->m_Points.rows(); // The number of vertices
  int numF = this->m_Faces.rows(); // The number of faces

  // Determine the faces containing each query point
  Eigen::VectorXi inFace( numQ );
  this->inElement( Xq, Yq, inFace );

  // Iterate over query points
  for( int i = 0; i < numQ; i++ ) {

    // The face containing the current query point
    int seedFace =  inFace(i);

    /*
    // Check whether the current query point lies within the convex hull of the
    // extended triangulation
    if ( seedFace == -1 ) {
      std::cout << "seedFace Check" << std::endl;
      std::runtime_error(
          "Query points lie outside convex hull of control points!");
    }
    */

    //-----------------------------------------------------------------------------------
    // Handle the case that a query point is equal to a control point
    //-----------------------------------------------------------------------------------
    
    // Find the shortest distance between the current query point and the control
    // points comprising the vertices of the containing face
    Scalar dx = this->m_Points(this->m_Faces(seedFace, 0), 0) - Xq(i);
    Scalar dy = this->m_Points(this->m_Faces(seedFace, 0), 1) - Yq(i);
    Scalar minDist = std::sqrt( dx*dx + dy*dy );
    int minID = this->m_Faces(seedFace, 0);

    for( int j = 1; j < 3; j++ ) {

      dx = this->m_Points(this->m_Faces(seedFace, j), 0) - Xq(i);
      dy = this->m_Points(this->m_Faces(seedFace, j), 1) - Yq(i);
      Scalar curDist = std::sqrt( dx*dx + dy*dy );

      if ( curDist < minDist ) {
        minDist = curDist;
        minID = this->m_Faces(seedFace, j);
      }

    }

    if ( minDist < Scalar(1.0e-12) ) {

      u[i] = Vector::Constant(1, 1, Scalar(1.0));
      uIDx[i] = Eigen::VectorXi::Constant(1, 1, minID);
      uVC[i] = Eigen::Matrix<Scalar, Eigen::Dynamic, 2>::Constant(1, 2, Scalar(0.0));
      uA(i) = Scalar(0.0);

      continue;

    }

    //-----------------------------------------------------------------------------------
    // Process faces
    //-----------------------------------------------------------------------------------
    // The union of the faces that contain the query point within their circumcircle form
    // a star-shaped 'natural neighbor polygon'. The vertices of this polygon are the
    // natural neighbors of the query point
    
    // A vector holding the area contributions to each control point
    // (Control points that are not natural neighbors of the current query point
    // will have zero contribution)
    Vector allQU = Vector::Constant( numP, 1, Scalar(0.0) );

    // A running sum to calculate the area of the virtual Voronoi polygon
    Scalar qA = 0.0;

    // A list of the CCW edges comporising the boundary of the
    // natural neighbor polygon
    std::vector<Eigen::RowVector2i> edgeList;

    // The real circumcenter associated to those edges
    std::vector<RowVec2d> realCC;

    // The virtual circumcenter associated to those edges
    std::vector<RowVec2d> virtualCC;

    // We determine the faces comprising the natural neighbor polygon by
    // performing a breadth-first search of the face graph structure starting
    // with the face containing the query point
    VectorXb visited = VectorXb::Constant( numF, 1, false );
    VectorXb inCircle = VectorXb::Constant( numF, 1, false );
    std::list<int> queue;

    visited(seedFace) = true;
    inCircle(seedFace) = true;
    queue.push_back(seedFace);

    while( !queue.empty() ) {

      // The ID of the current face
      int curF = queue.front();

      // The shifted circumcenter of the current face
      RowVec2d curFCC = this->m_FCC.row(curF);
      curFCC(0) = curFCC(0) - Xq(i);
      curFCC(1) = curFCC(1) - Yq(i);

      // Remove the current face from the queue
      queue.pop_front();

      // Iterate over the neighbors of the current face
      for( int j = 0; j < 3; j++ ) {

        int curFN = this->m_FaceNeighbors(curF, j);

        // Accumulate oriented boundary edges into the edge list
        if ( curFN == -1 ) {

          // Add the oriented edge to the edge list
          Eigen::RowVector2i curEdge = this->m_OrientedEdges.row(curF + numF*j);
          edgeList.push_back( curEdge );

          // Add the shifted circumcenter of the current face
          realCC.push_back( curFCC );

          // Calculate the shifted circumcenter of the triangle formed by
          // the query point and the two vertices defining the current edge
          RowVec2d curVCC = NNIpp::triCircumcenter( Xq(i), Yq(i),
              this->m_Points(curEdge(0), 0), this->m_Points(curEdge(0), 1),
              this->m_Points(curEdge(1), 0), this->m_Points(curEdge(1), 1) );

          curVCC(0) = curVCC(0) - Xq(i);
          curVCC(1) = curVCC(1) - Yq(i);

          virtualCC.push_back( curVCC );

          continue;

        }
        
        // Accumulate natural neighbor faces into the BFS queue
        if ( visited(curFN) == false ) {

          visited(curFN) = true;

          // Determine if the unvisited neighbor face is a member of
          // the natural neighbor polygon
          if ( this->inCircle(Xq(i), Yq(i), curFN) ) {
            inCircle(curFN) = true;
            queue.push_back( curFN );
          }

        }

        if ( inCircle(curFN) ) {

          // The shifted circumcenter of the current neighbor
          RowVec2d curFNCC = this->m_FCC.row(curFN);
          curFNCC(0) = curFNCC(0) - Xq(i);
          curFNCC(1) = curFNCC(1) - Yq(i);

          // Calculate the contribution to the control point area
          Scalar nodeArea = curFCC(0) * curFNCC(1) - curFCC(1) * curFNCC(0);

          allQU(this->m_Faces(curF, (j+2) % 3)) += nodeArea;

        } else {

          // Accumulate edges shared with non-natural neighbor faces into
          // the edge list
          Eigen::RowVector2i curEdge = this->m_OrientedEdges.row(curF + numF*j);
          edgeList.push_back( curEdge );

          // Add the shifted circumcenter of the current face
          realCC.push_back( curFCC );

          // Calculate the shifted circumcenter of the triangle formed by
          // the query point and the two vertices defining the current edge
          RowVec2d curVCC = NNIpp::triCircumcenter( Xq(i), Yq(i),
              this->m_Points(curEdge(0), 0), this->m_Points(curEdge(0), 1),
              this->m_Points(curEdge(1), 0), this->m_Points(curEdge(1), 1) );

          curVCC(0) = curVCC(0) - Xq(i);
          curVCC(1) = curVCC(1) - Yq(i);

          virtualCC.push_back( curVCC );

        }

      }

    }

    //-----------------------------------------------------------------------------------
    // Sort the edge list and associated lists of circumcenters to extract the natural
    // neighbors of the current query point and the vertices of its virtual Voronoi cell
    //-----------------------------------------------------------------------------------
    
    // The number of natural neighbors
    int numNN = edgeList.size();

    Eigen::Matrix<int, Eigen::Dynamic, 2> sortedEdges( numNN, 2 );
    Eigen::Matrix<Scalar, Eigen::Dynamic, 2> sortedRCC( numNN, 2 );
    Eigen::Matrix<Scalar, Eigen::Dynamic, 2> sortedVCC( numNN, 2 );

    sortedEdges.row(0) = edgeList.back();
    sortedRCC.row(0) = realCC.back();
    sortedVCC.row(0) = virtualCC.back();

    edgeList.pop_back();
    realCC.pop_back();
    virtualCC.pop_back();

    int count = 0;
    while( !edgeList.empty() ) {

      int vb = sortedEdges(count, 1);
      count++;

      for( int j = 0; j < edgeList.size(); j++ ) {

        Eigen::RowVector2i curEdge = edgeList[j];
        if ( curEdge(0) == vb ) {

          sortedEdges.row(count) = edgeList[j];
          sortedRCC.row(count) = realCC[j];
          sortedVCC.row(count) = virtualCC[j];

          std::swap( edgeList[j], edgeList.back() );
          std::swap( realCC[j], realCC.back() );
          std::swap( virtualCC[j], virtualCC.back() );

          edgeList.pop_back();
          realCC.pop_back();
          virtualCC.pop_back();

        }

      }

    }

    //-----------------------------------------------------------------------------------
    // Process area contributions to natural neighbor coordinates from each edge
    //-----------------------------------------------------------------------------------
    
    for( int j = 0; j < numNN; j++ ) {

      RowVec2d curVCC = sortedVCC.row(j);
      RowVec2d curRCC = sortedRCC.row(j);
      Scalar nodeArea = curVCC(0) * curRCC(1) - curVCC(1) * curRCC(0);

      allQU( sortedEdges(j,0) ) += nodeArea;
      allQU( sortedEdges(j,1) ) -= nodeArea;

      RowVec2d nextVCC = sortedVCC.row( (j+1) % numNN );
      nodeArea = curVCC(0) * nextVCC(1) - curVCC(1) * nextVCC(0);

      allQU( sortedEdges(j,1) ) += nodeArea;
      qA += nodeArea;

    }

    //-----------------------------------------------------------------------------------
    // Update output lists
    //-----------------------------------------------------------------------------------
    
    // The natural neighbors of the current query point
    Eigen::VectorXi quIDx = sortedEdges.col(0);

    // The natural neighbor coordinates of the current query point
    Vector qu(numNN);
    for( int j = 0; j < numNN; j++ ) {
      qu(j) = allQU(quIDx(j));
    }

    // Normalize coordinates
    qu = (qu.array() / qA).matrix();

    u[i] = qu;
    uIDx[i] = quIDx;
    uVC[i] = sortedVCC;
    uA(i) = qA / Scalar(2.0);

  }

};

// ======================================================================================
// INTERPOLATION FUNCTIONS
// ======================================================================================

///
/// Interpolate function at a set of query points using Sibson's
/// C^1 continuous method for scattered data interpolation
///
template <typename Scalar>
void NNIpp::NaturalNeighborInterpolant<Scalar>::operator()(
    const Vector &Xq, const Vector &Yq, Matrix &Fq ) const {

  typedef Eigen::Matrix<Scalar, 1, 2> RowVec2d;
  typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> ArrayVec;
  
  // The number of query points
  int numQ = Xq.size();

  // The number of function components
  int numVals = this->m_Values.cols();

  // Calculate the natural neighbor coordinates for each query point
  std::vector<Vector> allU( numQ );
  std::vector<Eigen::VectorXi> allUIDx( numQ );
  std::vector<Eigen::Matrix<Scalar, Eigen::Dynamic, 2> > allUVC( numQ );
  Vector allUA( numQ );

  this->naturalNeighborCoordinates( Xq, Yq, allU, allUIDx, allUVC, allUA );

  // Iterate over query points
  for( int i = 0; i < numQ; i++ ) {

    // The coordinates of the query point
    RowVec2d Vq;
    Vq << Xq(i), Yq(i);

    // The natural neighbor coordinates
    ArrayVec u = allU[i];

    // The IDs of the natural neighbors
    Eigen::ArrayXi uIDx = allUIDx[i];

    // The number of natural neighbors
    int numNN = u.size();

    // Handle the case where a query point coincides with a control point
    if ( numNN == 1 ) {

      for( int j = 0; j < numVals; j++ ) {
        Fq(i,j) = this->m_Values(uIDx(0), j);
      }

      continue;

    }

    // Calculate the separation vector/distances from the current query
    // point to each of its natural neighbors
    Eigen::Array<Scalar, Eigen::Dynamic, 2> d(numNN, 2);
    ArrayVec L(numNN, 1);

    for( int j = 0; j < numNN; j++ ) {

      d.row(j) = Vq - this->m_Points.row(uIDx(j));
      L(j) = std::sqrt(d(j,0)*d(j,0) + d(j,1)*d(j,1));
    
    }

    // Cache some convenience variables to avoid duplicate computations
    ArrayVec uLVec = u * L;
    ArrayVec uDLVec = u / L;
    Scalar uDLW = uDLVec.sum();

    // Calculate the 'alpha' and 'beta' parameters
    Scalar beta = ( uLVec * L ).sum();
    Scalar alpha = uLVec.sum() / uDLW;
    Scalar ab = alpha + beta;

    // Iterate over function components
    for( int j = 0; j < numVals; j++ ) {

      // Calculate the 'xi' parameters and the
      // C^0 interpolant
      ArrayVec xi_I(numNN, 1);
      Scalar Z0 = Scalar(0.0);

      for( int k = 0; k < numNN; k++ ) {

        xi_I(k) = this->m_Values(uIDx(k), j) +
          this->m_DataGradX(uIDx(k), j) * d(k,0) +
          this->m_DataGradY(uIDx(k), j) * d(k,1);

        Z0 += u(k) * this->m_Values(uIDx(k), j);

      }

      Scalar xi = ( uDLVec * xi_I ).sum() / uDLW;

      // Calculate the C^1 interpolant
      Fq(i,j) = ( alpha * Z0 + beta * xi ) / ab;

    }

  }

};

///
/// Interpolate function at a set of query points using Sibson's
/// C^1 continuous method for scattered data interpolation.
/// Also evaluate analytic function gradients
///
template <typename Scalar>
void NNIpp::NaturalNeighborInterpolant<Scalar>::operator()(
    const Vector &Xq, const Vector &Yq,
    Matrix &Fq, Matrix &DFx, Matrix &DFy ) const {

  typedef Eigen::Array<Scalar, 1, 2> RowVec2d;
  typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> ArrayVec;

  // The number of query points
  int numQ = Xq.size();

  // The number of function components
  int numVals = this->m_Values.cols();

  // Calculate the natural neighbor coordinates for each query point
  std::vector<Vector> allU( numQ );
  std::vector<Eigen::VectorXi> allUIDx( numQ );
  std::vector<Eigen::Matrix<Scalar, Eigen::Dynamic, 2> > allUVC( numQ );
  Vector allUA( numQ );

  this->naturalNeighborCoordinates( Xq, Yq, allU, allUIDx, allUVC, allUA );

  // Iterate over query points
  for( int i = 0; i < numQ; i++ ) {

    // The coordinates of the query point
    RowVec2d Vq;
    Vq << Xq(i), Yq(i);
    
    // The natural neighbor coordinates
    ArrayVec u = allU[i];

    // The IDs of the natural neighbors
    Eigen::ArrayXi uIDx = allUIDx[i];

    // The shifted vertices of the virtual Voronoi cell
    Eigen::Array<Scalar, Eigen::Dynamic, 2> uVC = allUVC[i];

    // The area of the virtual Voronoi cell
    Scalar uA = allUA(i);

    // The number of natural neighbors
    int numNN = u.size();

    // Handle the case where a query point conincides with a control point
    if ( numNN == 1 ) {

      for( int j = 0; j < numVals; j++ ) {
        Fq(i,j) = this->m_Values(uIDx(0), j);
        DFx(i,j) = this->m_DataGradX(uIDx(0), j);
        DFy(i,j) = this->m_DataGradY(uIDx(0), j);
      }

      continue;

    }

    // ----------------------------------------------------------------------------------
    // Calculate the gradients of the natural neighbor coordinates
    // ----------------------------------------------------------------------------------
    
    // The edge lengths of the virtual voronoi cell
    ArrayVec cellLengths( numNN, 1 );

    // The directed edges from natural neighbor control points to the
    // query points. NOTE: The gradient of the squared edge length is
    // simply the directed edge vector - hence the name
    Eigen::Array<Scalar, Eigen::Dynamic, 2> gradL2( numNN, 2 );

    // The squared lengths of those edges
    ArrayVec L2( numNN, 1 );

    // The lengths of those edges
    ArrayVec L( numNN, 1 );

    // The gradients of the edge lengths
    Eigen::Array<Scalar, Eigen::Dynamic, 2> gradL( numNN, 2 );

    // The gradients of the NON-NORMALIZED natural neighbor coordinates
    // with respect to the query point coordinates.
    Eigen::Array<Scalar, Eigen::Dynamic, 2> gradLU( numNN, 2 );

    int jPrev = numNN - 1;
    for( int j = 0; j < numNN; j++ ) {

      gradL2.row(j) = Vq - this->m_Points.row(uIDx(j)).array();

      L2(j) = gradL2(j,0)*gradL2(j,0) + gradL2(j,1)*gradL2(j,1);
      L(j) = std::sqrt(L2(j));

      gradL.row(j) << ( gradL2(j,0) / (Scalar(2.0)*L(j)) ),
        ( gradL2(j,1) / (Scalar(2.0)*L(j)) );

      RowVec2d cellVec = uVC.row(j) - uVC.row(jPrev);
      cellLengths(j) = std::sqrt( cellVec(0)*cellVec(0) + cellVec(1)*cellVec(1) );

      gradLU.row(j) = Scalar(0.5) * (uVC.row(j) + uVC.row(jPrev));
      gradLU.row(j) = cellLengths(j) * gradLU.row(j) / L(j);

      jPrev = j;

    }

    // The gradients of the NORMALIZED natural neighbor coordinates
    // with respect to the query point coordinates
    Eigen::Array<Scalar, Eigen::Dynamic, 2> gradU( numNN, 2 );
    gradU.col(0) = gradLU.col(0) - gradLU.col(0).sum() * u;
    gradU.col(1) = gradLU.col(1) - gradLU.col(1).sum() * u;
    gradU = gradU / uA;

    // ----------------------------------------------------------------------------------
    // Calculate the 'Beta' parameter and its derivatives
    // ----------------------------------------------------------------------------------
    
    Scalar beta = (u * L2).sum();

    RowVec2d gradBeta;
    gradBeta(0) = ( gradU.col(0) * L2 + u * gradL2.col(0) ).sum();
    gradBeta(1) = ( gradU.col(1) * L2 + u * gradL2.col(1) ).sum();

    // ----------------------------------------------------------------------------------
    // Calculate the 'Alpha' parameter and its derivatives
    // ----------------------------------------------------------------------------------
    
    ArrayVec w = u / L;
    Scalar WW = w.sum();

    Eigen::Array<Scalar, Eigen::Dynamic, 2> gradW(numNN, 2);
    gradW.col(0) = gradU.col(0) / L - u * gradL.col(0) / L2;
    gradW.col(1) = gradU.col(1) / L - u * gradL.col(1) / L2;

    RowVec2d gradWW;
    gradWW << (gradW.col(0).sum()), (gradW.col(1).sum());

    Scalar alpha = (u * L).sum() / WW;

    RowVec2d gradAlpha;
    gradAlpha(0) = ( L * gradU.col(0) + u * gradL.col(0) ).sum() - alpha * gradWW(0);
    gradAlpha(1) = ( L * gradU.col(1) + u * gradL.col(1) ).sum() - alpha * gradWW(1);
    gradAlpha = gradAlpha / WW;

    Scalar AB = alpha + beta;

    // Iterate over function components
    for( int j = 0; j < numVals; j++ ) {

      // --------------------------------------------------------------------------------
      // Calculate Sibson's C^0 interpolant and its derivatives
      // --------------------------------------------------------------------------------
      
      Scalar Z0 = Scalar(0.0);
      RowVec2d gradZ0 = RowVec2d::Zero();
      
      for( int k = 0; k < numNN; k++ ) {

        Z0 += u(k) * this->m_Values(uIDx(k), j);

        gradZ0(0) += gradU(k,0) * this->m_Values(uIDx(k), j);
        gradZ0(1) += gradU(k,1) * this->m_Values(uIDx(k), j);

      }

      // --------------------------------------------------------------------------------
      // Calculate the 'Xi' parameter and its derivatives
      // --------------------------------------------------------------------------------
      
      ArrayVec xi( numNN, 1 );
      Eigen::Array<Scalar, Eigen::Dynamic, 2> gradxi( numNN, 2 );
      for( int k = 0; k < numNN; k++ ) {

        xi(k) = this->m_Values(uIDx(k), j) +
          this->m_DataGradX(uIDx(k), j) * gradL2(k,0) +
          this->m_DataGradY(uIDx(k), j) * gradL2(k,1);

        gradxi(k,0) = this->m_DataGradX(uIDx(k), j);
        gradxi(k,1) = this->m_DataGradY(uIDx(k), j);

      }

      Scalar Xi = ( w * xi ).sum() / WW;

      RowVec2d gradXi;
      gradXi(0) = (gradW.col(0) * xi + w * gradxi.col(0)).sum() - Xi * gradW.col(0).sum();
      gradXi(1) = (gradW.col(1) * xi + w * gradxi.col(1)).sum() - Xi * gradW.col(1).sum();
      gradXi = gradXi / WW;

      // --------------------------------------------------------------------------------
      // Calculate Sibson's C^1 interpolant and its derivatives
      // --------------------------------------------------------------------------------
      
      Scalar Z1 = ( alpha * Z0 + beta * Xi ) / AB;

      RowVec2d gradZ1 = gradAlpha * (Z0-Z1) + gradBeta * (Xi-Z1) +
        alpha * gradZ0 + beta * gradXi;
      gradZ1 = gradZ1 / AB;

      // Update output variables
      Fq(i,j) = Z1;
      DFx(i,j) = gradZ1(0);
      DFy(i,j) = gradZ1(1);

    }

  }

};
      
//TODO: Add explicit template instantiation
#ifdef NNI_STATIC_LIBRARY
#endif
