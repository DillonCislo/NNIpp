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

#include "ghostPointsEdge.h"
#include "../General/delaunayTriangulation.h"

#include <cassert>
#include <stdexcept>
#include <iostream>

#include <igl/boundary_loop.h>

template <typename Scalar>
NNI_INLINE void NNIpp::ghostPointsEdge(
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &Xp,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &Yp,
    const NNIParam<Scalar> &param,
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &GPx,
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &GPy ) {

  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
  typedef Eigen::Matrix<Scalar, 2, 1> Vector2d;

  // Construct the Delaunay triangulation of the input points
  Eigen::Matrix<int, Eigen::Dynamic, 3> F;
  NNIpp::delaunayTriangulation( Xp, Yp, F );

  // Determine the boundary vertices of the triangulation
  Eigen::VectorXi bdyIDx;
  igl::boundary_loop( F, bdyIDx );

  int numBdyV = bdyIDx.rows();
  GPx = Vector::Zero( numBdyV );
  GPy = Vector::Zero( numBdyV );

  for( int i = 0; i < numBdyV; i++ ) {

    int j = (i+1) % numBdyV;

    int v1 = bdyIDx(i);
    int v2 = bdyIDx(j);

    // Calculate the edge midpoint
    Scalar EMx = (Xp(v1)+Xp(v2)) / Scalar(2.0);
    Scalar EMy = (Yp(v1)+Yp(v2)) / Scalar(2.0);

    // Calculate the CW rotated edge vectors
    Scalar EVx = param.GPe * (Yp(v2)-Yp(v1));
    Scalar EVy = param.GPe * (Xp(v1)-Xp(v2));

    // Calculate the ghost point location
    GPx(i) = EMx + EVx;
    GPy(i) = EMy + EVy;

  }

};

// TODO: Add explicit template instantiation
#ifdef NNI_STATIC_LIBRARY
#endif
