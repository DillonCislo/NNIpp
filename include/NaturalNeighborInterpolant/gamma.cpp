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

#include "gamma.h"

#include <cassert>
#include <cmath>
#include <iostream>

template <typename Scalar>
void NNIpp::gamma(
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 4> &X,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 4> &Y,
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &G ) {

  typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> Vector;

  // Check that the input arguments are properly sized
  assert( X.rows() == Y.rows() );
  assert( G.rows() == X.rows() );

  // Extract the x-coordinates of each vertex
  Vector x1 = X.col(0).array();
  Vector x2 = X.col(1).array();
  Vector x3 = X.col(2).array();
  Vector x4 = X.col(3).array();

  // Extract the y-coordinates of each vertex
  Vector y1 = Y.col(0).array();
  Vector y2 = Y.col(1).array();
  Vector y3 = Y.col(2).array();
  Vector y4 = Y.col(3).array();

  // Cache some variables to avoid duplicate computations
  Vector x1S = x1 * x1;
  Vector x2S = x2 * x2;
  Vector x3S = x3 * x3;
  Vector x4S = x4 * x4;

  Vector y1S = y1 * y1;
  Vector y2S = y2 * y2;
  Vector y3S = y3 * y3;
  Vector y4S = y4 * y4;

  Vector x1y4y3 = x1 * (y4 - y3);
  Vector x1y2y4 = x1 * (y2 - y4);
  Vector x1y3y2 = x1 * (y3 - y2);

  Vector x2y3y4 = x2 * (y3 - y4);
  Vector x2y4y1 = x2 * (y4 - y1);
  Vector x2y1y3 = x2 * (y1 - y3);

  Vector x3y4y2 = x3 * (y4 - y2);
  Vector x3y1y4 = x3 * (y1 - y4);
  Vector x3y2y1 = x3 * (y2 - y1);

  Vector x4y2y3 = x4 * (y2 - y3);
  Vector x4y3y1 = x4 * (y3 - y1);
  Vector x4y1y2 = x4 * (y1 - y2);

  Vector c1S = x2y3y4 + x3y4y2 + x4y2y3;
  Vector c2S = x1y4y3 + x3y1y4 + x4y3y1;
  Vector c3S = x1y2y4 + x2y4y1 + x4y1y2;
  Vector c4S = x1y3y2 + x2y1y3 + x3y2y1;

  // Construct the 'Gamma' vector
  Vector GArr = c1S * ( x1S + y1S ) + c2S * ( x2S + y2S ) +
                c3S * ( x3S + y3S ) + c4S * ( x4S + y4S );

  G = GArr.matrix();

};
