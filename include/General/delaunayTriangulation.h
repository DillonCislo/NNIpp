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

#ifndef _DELAUNAY_TRIANGULATION_H_
#define _DELAUNAY_TRIANGULATION_H_

#include <Eigen/Core>

namespace NNIpp {

  ///
  /// Construct the Delaunay triangulation of a set of input points
  /// This function has an unfortunate dependency on CGAL. This will
  /// be removed once the libigl Delaunay triangulation is fixed and
  /// made more robust
  ///
  /// Templates:
  ///
  ///   Scalar    Input type of scattered data points, values, and gradients
  ///
  /// Inputs:
  ///
  ///   X     #P by 1 list of x-coordinates
  ///   Y     #P by 1 list of y-coordinates
  ///
  /// Outputs:
  ///
  ///   F     #F by 3 face connectivity list
  ///
  template <typename Scalar>
  void delaunayTriangulation(
      const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &X,
      const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &Y,
      Eigen::Matrix<int, Eigen::Dynamic, 3> &F );

}

#endif
