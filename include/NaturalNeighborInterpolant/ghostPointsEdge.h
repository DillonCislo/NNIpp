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

#ifndef _GHOST_POINTS_EDGE_H_
#define _GHOST_POINTS_EDGE_H_

#include "../General/nniInline.h"
#include <Eigen/Core>

namespace NNIpp {

  ///
  /// Generate the ghost points used for natural neighbor extrapolation by
  /// using the outward rotated mid-edge vectors of the Delaunay triangulation
  /// of the input scattered data points
  ///
  /// Templates:
  ///
  ///   Scalar    Input type of scattered data points, values, and gradients
  ///
  /// Inputs:
  ///
  ///   Xp        #P by 1 list of data point x-coordinates
  ///   Yp        #P by 1 list of data point y-coordinates
  ///   param     An 'NNIParam' class containing the rest of the parameters
  ///             needed to construct the interpolant
  ///
  /// Outputs:
  ///
  ///   GPx       #GP by 1 list of ghost point x-coordinates
  ///   GPy       #GP by 1 list of ghost point y-coordinates
  ///
  template <typename Scalar>
  NNI_INLINE void ghostPointsEdge(
      const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &Xp,
      const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &Yp,
      const NNIParam<Scalar> &param,
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &GPx,
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &GPy );

}

#ifndef NNI_STATIC_LIBRARY
#  include "ghostPointsEdge.cpp"
#endif

#endif
