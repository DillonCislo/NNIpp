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

#ifndef _GRAPH_DISTANCES_H_
#define _GRAPH_DISTANCES_H_

#include "nniInline.h"
#include <Eigen/Core>

namespace NNIpp {

  ///
  /// Calculates the length of the shortest paths between all pairs of
  /// vertices of a mesh triangulation along mesh edges using a modified
  /// Dijkstra's algorithm
  ///
  /// Templates:
  ///
  ///   DerivedV    Container type for mesh vertices
  ///   DerivedF    Container type for mesh connectivity list
  ///   DerivedD    Container type for graph distance matrix
  ///
  /// Inputs:
  ///
  ///   V   #V by dim vertex coordinate list
  ///   F   #F by 3 face connectivity list
  ///
  /// Outputs:
  ///
  ///   D   #V by #V matrix. D(i,j) is the length of the shortest path
  ///       between node i and node j
  ///
  template <typename DerivedV, typename DerivedF, typename DerivedD>
  NNI_INLINE void graphDistances(
      const Eigen::MatrixBase<DerivedV> &V,
      const Eigen::MatrixBase<DerivedF> &F,
      Eigen::MatrixBase<DerivedD> &D );

  ///
  /// Calculates the length of the shortest paths between all pairs of
  /// vertices of a mesh triangulation along unit weighted mesh edges
  /// using a modified breadth-first search
  ///
  /// Templates:
  ///
  ///   DerivedF    Container type for mesh connectivity list
  ///   DerivedD    Container type for graph distance matrix
  ///
  /// Inputs:
  ///
  ///   F   #F by 3 face connectivity list
  ///
  /// Outputs:
  ///
  ///   D   #V by #V matrix. D(i,j) is the length of the shortest path
  ///       between node i and node j
  ///
  template <typename DerivedF, typename DerivedD>
  NNI_INLINE void graphDistances(
      const Eigen::MatrixBase<DerivedF> &F,
      Eigen::MatrixBase<DerivedD> &D );

}

#ifndef NNI_STATIC_LIBRARY
#  include "graphDistances.cpp"
#endif

#endif
