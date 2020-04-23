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

#ifndef _GRADIENTS_DIRECT_H_
#define _GRADIENTS_DIRECT_H_

#include "../General/nniInline.h"
#include "NNIParam.h"
#include <Eigen/Core>

namespace NNIpp {

  ///
  /// Generate discrete gradient/Hessian data for scattered input data
  /// Discrete derivatives are the 1st- and 2nd-order terms of a 3rd-
  /// order Talyor polynomial fit to the 3rd-order natural neighborhood
  /// of each input point
  ///
  /// Templates:
  ///
  ///   Scalar    Input type of scattered data points, values, and gradients
  ///
  /// Inputs:
  ///
  ///   V       (#P+#GP) by 2 list of coordinates of the extended triangulation
  ///   F       #F by 3 face connectivity list of the extended triangulation
  ///   Vals    (#P+#GP) by #V matrix of values for the extended triangulation
  ///   param   NNIParam class. Contains parameters for interpolant construction
  ///
  /// Outputs:
  ///
  ///   Dx      (#P+#GP) by #V matrix of x-gradient values
  ///   Dy      (#P+#GP) by #V matrix of y-gradient values
  ///   Dxx     (#P+#GP) by #V matrix of (x,x)-Hessian entries
  ///   Dxy     (#P+#GP) by #V matrix of (x,y)-Hessian entries
  ///   Dyy     (#P+#GP) by #V matrix of (y,y)-Hessian entries
  ///
  template <typename Scalar>
  NNI_INLINE void gradientsDirect(
      const Eigen::Matrix<Scalar, Eigen::Dynamic, 2> &V,
      const Eigen::Matrix<int, Eigen::Dynamic, 3> &F,
      const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &Vals,
      const NNIParam<Scalar> &param,
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &Dx,
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &Dy,
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &Dxx,
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &Dxy,
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &Dyy );

}

#ifndef NNI_STATIC_LIBRARY
#  include "gradientsDirect.cpp"
#endif

#endif
