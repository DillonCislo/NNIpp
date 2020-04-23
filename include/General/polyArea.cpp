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

#include "polyArea.h"

template <typename Scalar>
NNI_INLINE Scalar NNIpp::polyArea(
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &Px,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &Py ) {

  // The number of polygon vertices
  int numV = Px.size();

  // The polygon area
  Scalar A = Scalar(0.0);

  int j = numV-1;
  for( int i = 0; i < numV; i++ ) {

    A += ( Px(i) + Px(j) ) * ( Py(i) - Py(j) );
    j = i;

  }

  return ( A / Scalar(2.0) );

};

template <typename Scalar>
NNI_INLINE Scalar NNIpp::polyArea(
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 2> &P ) {

  // The number of polygon vertices
  int numV = P.rows();

  // The polygon area
  Scalar A = Scalar(0.0);

  int j = numV-1;
  for( int i = 0; i < numV; i++ ) {

    A += ( P(i,0) + P(j,0) ) * ( P(i,1) - P(j,1) );
    j = i;

  }

  return ( A / Scalar(2.0) );

};

//TODO: Add explicit template instantiation
#ifdef NNI_STATIC_LIBRARY
#endif


 
