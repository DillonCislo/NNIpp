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

#ifndef _IN_POLYGON_H_
#define _IN_POLYGON_H_

#include "nniInline.h"
#include <Eigen/Core>

namespace NNIpp {

  ///
  /// Determine if a query point lies within a planar polygon
  /// WARNING: This method is fast, but is unreliable for points
  /// that lie directly on polygon edges
  ///
  /// Templates:
  ///
  ///   Scalar    Input type of scattered data points, values, and gradients
  ///
  /// Inputs:
  ///
  ///   poly    #P by 2 ordered list of polygon vertex (x,y) coordinates
  ///   Xq      #Q by 1 list of query point x-coordinates
  ///   Yq      #Q by 1 list of query point y-coordinates
  ///
  /// Outputs:
  ///
  ///   inPoly  #Q by 1 boolean array. Entries are true if the
  ///           corresponding query point lies within or on the
  ///           boundary of the input polygon
  ///
  template <typename Scalar>
  NNI_INLINE void inPolygon(
      const Eigen::Matrix<Scalar, Eigen::Dynamic, 2> &poly,
      const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &Xq,
      const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &Yq,
      Eigen::Array<bool, Eigen::Dynamic, 1> &inPoly );

};

#ifndef NNI_STATIC_LIBRARY
#  include "inPolygon.cpp"
#endif

#endif
