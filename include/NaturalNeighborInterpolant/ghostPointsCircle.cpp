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

#include "ghostPointsCircle.h"

#include <cassert>
#include <cmath>
#include <iostream>

template <typename Scalar>
NNI_INLINE void NNIpp::ghostPointsCircle(
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &Xp,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &Yp,
    const NNIpp::NNIParam<Scalar> &param,
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &GPx,
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &GPy ) {

  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
  typedef Eigen::Matrix<Scalar, 2, 1> Vector2d;

  // The bounding box of the input data
  Vector2d BBx;
  BBx << Xp.minCoeff(), Xp.maxCoeff();

  Vector2d BBy;
  BBy << Yp.minCoeff(), Yp.maxCoeff();

  // The side lengths of the bounding box
  Scalar LBBx = BBx(1)-BBx(0);
  Scalar LBBy = BBy(1)-BBy(0);

  // The center of the circumcircle of the bounding box 
  Scalar BBCx = BBx(0) + LBBx / Scalar(2.0);
  Scalar BBCy = BBy(0) + LBBy / Scalar(2.0);

  // The radius of the circumcircle of the bounding box
  Scalar BBCr;
  BBCr = std::sqrt( LBBx * LBBx + LBBy * LBBy ) / Scalar(2.0);
  BBCr = param.GPr * BBCr; // Increase radius

  // Calculate the locations of the ghost points
  Scalar size = (Scalar) (2.0 * M_PI / param.GPn);
  Scalar low = Scalar(0.0);
  Scalar high = (Scalar) ((2.0 * M_PI) - size);

  Vector theta(param.GPn);
  for( int i = 0; i < param.GPn; i++ ) {
    theta(i) = low + i * size;
  }

  GPx = theta.array().cos().matrix();
  GPx = ( (BBCr * GPx.array()) + BBCx ).matrix();

  GPy = theta.array().sin().matrix();
  GPy = ( (BBCr * GPy.array()) + BBCy ).matrix();

};

// TODO: Add explicit template instantiation
#ifdef NNI_STATIC_LIBRARY
#endif
