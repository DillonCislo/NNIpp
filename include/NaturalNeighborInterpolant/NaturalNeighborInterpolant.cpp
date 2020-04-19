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

#include <igl/edges.h>

#include "ghostPointsCircle.h"
#include "ghostPointsEdge.h"
#include "../General/delaunayTriangulation.h"
#include "../General/convexHull.h"
#include "../General/inPolygon.h"
#include "../General/triCircumcenter.h"
#include "delta.h"

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
  this->m_Points =
    Eigen::Matrix<Scalar, Eigen::Dynamic, 2>::Zero( (Xp.size() + GPx.size()), 2 );
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

  // Calculate the convex hull of the extended triangulation
  this->m_ConvexHull = Eigen::Matrix<Scalar, Eigen::Dynamic, 2>::Zero(1,2);
  NNIpp::convexHull( Px, Py, this->m_ConvexHull );

  // Calculate the 'delta' parameters for each face of the extended triangulation
  NNIpp::delta( FVx, FVy, this->m_fDelta );

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

};

//TODO: Add explicit template instantiation
#ifdef NNI_STATIC_LIBRARY
#endif
