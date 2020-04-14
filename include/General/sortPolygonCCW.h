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

#ifndef _SORT_POLYGON_CCW_H_
#define _SORT_POLYGON_CCW_H_

#include "nniInline.h"
#include <Eigen/Core>

namespace NNIpp {

  /// Sort the vertices of a convex planar polygon in counter-clockwise order
  ///
  /// Templates:
  ///
  ///   Scalar    Input type of scattered data points, values, and gradients
  ///
  /// Inputs:
  ///
  ///   poly      #P by 2 list of unordered polygon vertices
  ///
  /// Outputs:
  ///
  ///   ccwPoly   #P by 2 list of sorted polygon vertices
  ///   order     #P by 1 list of indices into the old ordering
  ///             that produces the new ordering
  ///
  template <typename Scalar>
  NNI_INLINE void sortPolygonCCW(
      const Eigen::Matrix<Scalar, Eigen::Dynamic, 2> &poly,
      Eigen::Matrix<Scalar, Eigen::Dynamic, 2> &ccwPoly,
      Eigen::VectorXi &order );

}

#ifndef NNI_STATIC_LIBRARY
#  include "sortPolygonCCW.cpp"
#endif

#endif
