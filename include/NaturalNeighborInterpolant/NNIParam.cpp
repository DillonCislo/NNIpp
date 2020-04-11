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

#include "NNIParam.h"

#include <cassert>
#include <stdexcept>
#include <iostream>

//
// Constructor for interpolation parameter class
//
NNIpp::NNIParam::NNIParam(int nPts) : numPoints(nPts) {

  ghostMethod = NNI_GHOST_POINTS_EDGE;
  GPe = Scalar(1.0);
  GPr = Scalar(2.0);
  GPn = 500;
  gradType = NNI_GRADIENTS_DIRECT;
  iterAlpha = Scalar(0.001);
  dispGrad = true;
  dispInterp = true;
  gradientSupplied = false;
  gradientSupplied = false;

};

//
// Check the validity of the interpolation parameters
//
inline void NNIpp::NNIParam::checkParam() const {

  if (numPoints <= 0)
    throw std::invalid_argument("Number of points must be positive");
  if ( ghostMethod < NNI_GHOST_POINTS_CUSTOM ||
      ghostMethod > NNI_GHOST_POINTS_EDGE )
    throw std::invalid_argument("Unsupported ghost point construction method");
  if ( GPe < Scalar(1.0) )
    throw std::invalid_argument("'GPe' must not be less than one");
  if ( GPr < Scalar(1.0) )
    throw std::invalid_argument("'GPr' must not be less than one");
  if ( GPn <= 0 )
    throw std::invalid_argument("'GPn' must be positive");
  if ( gradType < NNI_GRADIENTS_DIRECT ||
      gradType > NNI_GRADIENTS_ITER )
    throw std::invalid_argument("Unsupported gradient generation method");
  if ( iterAlpha < Scalar(0.0) || iterAlpha > Scalar(1.0) )
    throw std::invalid_argument("'interAlpha' must lie in the range [0,1]");

};
