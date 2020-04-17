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

#include "convexHull.h"

#include <cassert>
#include <stdexcept>
#include <iostream>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>

template <typename Scalar>
NNI_INLINE void NNIpp::convexHull(
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &X,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &Y,
    Eigen::Matrix<Scalar, Eigen::Dynamic, 2> &CH ) {

  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef Kernel::Point_2 Point_2;

  // Convert Eigen-style input coordinates to C++ style vector
  int numPoints = X.size();
  std::vector<Point_2> points( X.size() );
  for( int i = 0; i < numPoints; i++ ) {
    points[i] = Point_2( X(i), Y(i) );
  }

  // Determine the convex hull
  std::vector<Point_2> result;
  CGAL::convex_hull_2( points.begin(), points.end(),
      std::back_inserter(result) );

  // Format output
  CH.resize( result.size(), 2 );
  for( int i = 0; i < result.size(); i++ ) {
    CH(i,0) = result[i].x();
    CH(i,1) = result[i].y();
  }

};

//TODO: Add explicit template instantiation
#ifdef NNI_STATIC_LIBRARY
#endif
