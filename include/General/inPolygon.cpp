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


#include "inPolygon.h"

#include <cassert>
#include <stdexcept>
#include <iostream>
#include <limits>

template <typename Scalar>
NNI_INLINE void NNIpp::inPolygon(
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 2> &poly,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &Xq,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &Yq,
    Eigen::Array<bool, Eigen::Dynamic, 1> &inPoly ) {

  typedef Eigen::Array<bool, Eigen::Dynamic, 1> ArrayXb;

  // The number of query points
  int numPoints = Xq.size();
  inPoly = ArrayXb::Constant( numPoints, false );

  // The number of vertices comprising the input polygon
  int numV = poly.rows();
  if (numV < 3)
    return;

  Scalar xold = poly(numV-1, 0);
  Scalar yold = poly(numV-1, 1);
  Scalar xnew, ynew;
  Scalar x1, y1;
  Scalar x2, y2;
  bool cX, cY;

  for( int i = 0; i < numV; i++ ) {

    // The coordinates of the current vertex
    xnew = poly(i,0);
    ynew = poly(i,1);

    // Construct the oriented edge vector
    if (xnew > xold) {

      x1 = xold;
      y1 = yold;
      x2 = xnew;
      y2 = ynew;

    } else {

      x1 = xnew;
      y1 = ynew;
      x2 = xold;
      y2 = yold;

    }

    // Update the edge crossing count for each query point
    for( int j = 0; j < numPoints; j++ ) {

      cX = (x1 < Xq(j)) && (Xq(j) <= x2);
      cY = ( (Yq(j)-y1) * (x2-x1) ) < ( (y2-y1) * (Xq(j)-x1) );

      if (cX && cY)
        inPoly(j) = !inPoly(j);

    }

    // Update the previous vertex
    xold = xnew;
    yold = ynew;

  }

};

// TODO: Add explicit template instantiation
#ifdef NNI_STATIC_LIBRARY
#endif
