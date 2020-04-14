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

#include "delaunayTriangulation.h"

#include <cassert>
#include <stdexcept>
#include <iostream>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

template <typename Scalar>
void NNIpp::delaunayTriangulation(
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &X,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &Y,
    Eigen::Matrix<int, Eigen::Dynamic, 3> &F ) {

  typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;
  typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, Kernel> Vb;
  typedef CGAL::Triangulation_data_structure_2<Vb>                Tds;
  typedef CGAL::Delaunay_triangulation_2<Kernel, Tds>             Delaunay;
  typedef Kernel::Point_2                                         Point;

  // The number of input points
  int numPoints = X.rows();

  // Convert Eigen-style input vectors to a C++ style vectors
  std::vector< std::pair<Point, unsigned> > points( numPoints );
  for( int i = 0; i < numPoints; i++ ) {

    points[i] = std::make_pair( Point(X(i), Y(i)), i );

  }

  // Construct Delaunay triangulation
  Delaunay triangulation;
  triangulation.insert( points.begin(), points.end() );

  // Format output as Eigen-style matrix
  F.resize( triangulation.number_of_faces(), 3 );

  int j = 0;
  for( Delaunay::Finite_faces_iterator fit = triangulation.finite_faces_begin();
      fit != triangulation.finite_faces_end(); ++fit ) {

    Delaunay::Face_handle face = fit;
    F(j,0) = face->vertex(0)->info();
    F(j,1) = face->vertex(1)->info();
    F(j,2) = face->vertex(2)->info();

    j++;

  }

};

// TODO: Add explicit template instantiation
#ifdef NNI_STATIC_LIBRARY
#endif
