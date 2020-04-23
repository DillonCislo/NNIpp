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

#ifndef _POLY_AREA_H_
#define _POLY_AREA_H_

#include "nniInline.h"
#include <Eigen/Core>

namespace NNIpp {

  ///
  /// Calculate the area of simple planar (convex) polygon
  /// NOTE: Polygon vertices must be counter-clockwise ordered
  ///
  /// Templates:
  ///
  ///   Scalar    Input type of scattered data points, values, and gradients
  ///
  /// Inputs:
  ///
  ///   Px      #P by 1 CCW sorted list of polygon vertex x-coordinates
  ///   Py      #P by 1 CCW sorted list of polygon vertex y-coordinates
  ///
  /// Outputs:
  ///
  ///   A       The enclosed area of the polygon
  ///
  template <typename Scalar>
  NNI_INLINE Scalar polyArea(
      const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &Px,
      const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &Py );

  ///
  /// Calculate the area of simple planar (convex) polygon
  /// NOTE: Polygon vertices must be counter-clockwise ordered
  ///
  /// Templates:
  ///
  ///   Scalar    Input type of scattered data points, values, and gradients
  ///
  /// Inputs:
  ///
  ///   P       #P by 2 CCW sorted list of polygon vertex (x,y)-coordinates
  ///
  /// Outputs:
  ///
  ///   A       The enclosed area of the polygon
  ///
  template <typename Scalar>
  NNI_INLINE Scalar polyArea(
      const Eigen::Matrix<Scalar, Eigen::Dynamic, 2> &P );

}

#ifndef NNI_STATIC_LIBRARY
#  include "polyArea.cpp"
#endif

#endif
