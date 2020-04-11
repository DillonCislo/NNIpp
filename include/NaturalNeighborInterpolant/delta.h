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

#ifndef _DELTA_H_
#define _DELTA_H_

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
  void delta(
      const Eigen::Matrix<Scalar, Eigen::Dynamic, 3> &X,
      const Eigen::Matrix<Scalar, Eigen::Dynamic, 3> &Y,
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &D );

}

#endif

