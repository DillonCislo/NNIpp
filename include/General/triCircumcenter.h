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

#ifndef _TRICIRCUMCENTER_H_
#define _TRICIRCUMCENTER_H_

#include "nniInline.h"
#include <Eigen/Core>

namespace NNIpp {

  ///
  /// A vectorized calculation of the circumcenters of a set of
  /// planar triangles
  ///
  /// Templates:
  ///
  ///   Scalar    Input type of scattered data points, values, and gradients
  ///
  /// Inputs:
  ///
  ///   X     #T by 3 list of triangle vertex x-coordinates
  ///   Y     #T by 3 list of triangle vertex y-coordinates
  ///
  /// Outputs:
  ///
  ///   CC    #T by 2 list of circumcenter coordinates
  ///
  template <typename Scalar>
  NNI_INLINE void triCircumcenter(
      const Eigen::Matrix<Scalar, Eigen::Dynamic, 3> &X,
      const Eigen::Matrix<Scalar, Eigen::Dynamic, 3> &Y,
      Eigen::Matrix<Scalar, Eigen::Dynamic, 2> &CC );

}

#ifndef NNI_STATIC_LIBRARY
#  include "triCircumcenter.cpp"
#endif

#endif
