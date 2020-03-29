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

#ifndef _CALC_DELTA_H_
#define _CALC_DELTA_H_

#include <cassert>
#include <stdexcept>
#include <iostream>

#include <Eigen/Core>

namespace NNIpp {

  ///
  /// Calculate the 'Delta' parameter from (Hiyoshi, 2008).
  /// Delta(v1, v2, v3) is equal to twice the signed area of the
  /// triangle defind by the points {v1, v2, v3}
  ///
  /// Templates:
  ///
  ///   Scalar    Input type of scattered data point, values, and gradients
  ///
  /// Inputs:
  ///
  ///   X         #N by 3 list of x-coordinates
  ///   Y         #N by 3 list of y-coordinates
  ///
  /// Outputs:
  ///
  ///   D         #N by 1 list of 'Delta' values
  ///
  template <typename Scalar>
  void calcDelta(
      const Eigen::Matrix<Scalar, Eigen::Dynamic, 3> &X,
      const Eigen::Matrix<Scalar, Eigen::Dynamic, 3> &Y,
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &D ) {

    typename Eigen::Array<Scalar, Eigen::Dynamic, 1> Vector;

    // Check that the input arguments are properly sized
    assert( X.rows() == Y.rows() );
    assert( X.rows() == D.rows() );

    // Extract the x-coordinates of each vertex
    Vector x1 = X.col(0).array();
    Vector x2 = X.col(1).array();
    Vector x3 = X.col(2).array();

    // Extract the y-coordinates of each vertex
    Vector y1 = Y.col(0).array();
    Vector y2 = Y.col(1).array();
    Vector y3 = Y.col(2).array();

    // Calculate 'Delta' parameter
    Vector Darr = x3 * ( y1 - y2 ) + x1 * ( y2 - y3 ) + x2 * ( y3 - y1 );

    D = Darr.matrix();

  };

} // namespace NNIpp

#endif // _CALC_DELTA_H_

