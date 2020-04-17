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

#ifndef _CONVEX_HULL_H_
#define _CONVEX_HULL_H_

#include "nniInline.h"
#include <Eigen/Core>

namespace NNIpp {

  ///
  /// Determine the convex hull of a set of points in the plane
  /// This function has an unfortunate dependency on CGAL. This will
  /// be removed once either the libigl developers or I choose to write
  /// one from scratch
  ///
  /// Templates:
  ///
  ///   Scalar    Input type of scattered data points, values, and gradients
  ///
  /// Inputs:
  ///
  ///   X   #P by 1 list of data point x-coordinates
  ///   Y   #P by 1 list of data point y-coordinates
  ///
  /// Outputs:
  ///
  ///   CH    #CH by 2 list of convex hull (x,y)-coordinates
  ///
  template <typename Scalar>
  NNI_INLINE void convexHull(
      const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &X,
      const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &Y,
      Eigen::Matrix<Scalar, Eigen::Dynamic, 2> &CH );

};

#ifndef NNI_STATIC_LIBRARY
#  include "convexHull.cpp"
#endif

#endif
