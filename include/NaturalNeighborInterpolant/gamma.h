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

#ifndef _GAMMA_H_
#define _GAMMA_H_

#include <Eigen/Core>

namespace NNIpp {

  ///
  /// Calculate the 'Gamma' parameter from (Hiyoshi, 2008). Used to
  /// calculate the natural neighbor coordinates of a query point
  /// Gamma(v1, v2, v3, v4) does not have a simple geometric interpretation
  ///
  /// Templates:
  ///
  ///   Scalar    Input type of scattered data points, values, and gradients
  ///
  /// Inputs:
  ///
  ///   X         #N by 4 list of x-coordinates
  ///   Y         #N by 4 list of y-coordinates
  ///
  /// Outputs:
  ///
  ///   G         #N by 1 list of 'Gamma' values
  ///
  template <typename Scalar>
  void gamma(
      const Eigen::Matrix<Scalar, Eigen::Dynamic, 4> &X,
      const Eigen::Matrix<Scalar, Eigen::Dynamic, 4> &Y,
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &G );

}

#endif
