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

#include <cassert>
#include <cmath>
#include <iostream>

#include <igl/sort_angles.h>

template <Scalar>
void NNIpp::sortPolygonCCW(
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 2> &poly,
    Eigen::Matrix<Scalar, Eigen::Dynamic, 2> &ccwPoly,
    Eigen::VectorXi &order ) {

  // The number of vertices in the polygon
  const std::size_t numRows = poly.rows();

  // Calculate the center of mass of the vertex coordinates
  Eigen::Matrix<Scalar, 1, 2> COMRow;
  COMRow << ( poly.col(0).mean() ), ( poly.col(1).mean() );

  Eigen::Matrix<Scalar, Eigen::Dynamic, 2> COM = COMRow.replicate( numRows, 1 );

  // Shift the polygon vertices so that the center of mass is at the origin
  Eigen::Matrix<Scalar, Eigen::Dynamic, 2> polyShift = poly - COM;

  // Sort the angles of the shifted vertices in ascending order
  igl::sort_angles( polyShift, order );

  // Re-order the polygon vertices
  for ( int i = 0; i < numRows; i++ ) {
    ccwPoly.row(i) = poly.row( order(i) );
  }

};
