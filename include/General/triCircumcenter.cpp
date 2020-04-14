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

#include "triCircumcenter.h"

#include <cassert>
#include <stdexcept>
#include <iostream>


template <typename Scalar>
NNI_INLINE void NNIpp::triCircumcenter(
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 3> &X,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 3> &Y,
    Eigen::Matrix<Scalar, Eigen::Dynamic, 2> &CC ) {

  typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> Vector;

  // Check that the input arguments are properly sized
  assert( (X.rows() == Y.rows()) && (X.cols() == Y.cols()) );
  assert( (X.rows() == CC.rows()) && (X.cols() == CC.cols()) );

  // Extract the x-coordinates of each vertex
  Vector x1 = X.col(0).array();
  Vector x2 = X.col(1).array();
  Vector x3 = X.col(2).array();

  // Extract the y-coordinates of each vertex
  Vector y1 = Y.col(0).array();
  Vector y2 = Y.col(1).array();
  Vector y3 = Y.col(2).array();

  // Cache some variables to avoid duplicate computations
  Vector y12 = y2-y1;
  Vector y31 = y1-y3;
  Vector y23 = y3-y2;

  Vector y12_P = y2+y1;
  Vector y31_P = y1+y3;
  Vector y23_P = y3+y2;

  Vector y123 = y12 * y23 * y31;

  Vector x1S = x1 * x1;
  Vector x2S = x2 * x2;
  Vector x3S = x3 * x3;

  Vector x12S = x2S - x1S;
  Vector x23S = x3S - x2S;
  Vector x31S = x1S - x3S;

  Vector x3y12 = x3 * y12;
  Vector x2y31 = x2 * y31;
  Vector x1y23 = x1 * y23;

  // The numerator for the x-coordinates
  Vector xNum = x3 * x3y12 + x2 * x2y31 + x1 * x1y23 - y123;

  // The numerator for the y-coordinates
  Vector yNum = x3 * x12S + x2 * x31S + x1 * x23S +
                x3y12 * y12_P + x2y31 * y31_P + x1y23 * y23_P;

  // The denominator for both coordinates
  Vector denom = Scalar(2.0) * ( x3y12 + x2y31 + x1y23 );

  // Calculate circumcenter coordinates
  Eigen::Array<Scalar, Eigen::Dynamic, 2> CCArr( X.rows(), 2 );
  CCArr << (xNum / denom), (yNum / denom);

  // Format output
  CC = CCArr.matrix();

};

// TODO: Add explicit template instantiation
#ifdef NNI_STATIC_LIBRARY
#endif
