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

#ifndef _GRAPH_NEIGHBORHOOD_H_
#define _GRAPH_NEIGHBORHOOD_H_

#include "nniInline.h"
#include <Eigen/Core>
#include <vector>

namespace NNIpp {

  ///
  /// Determine the N-th order neighborhood for each vertex in a mesh
  /// triangulation, i.e. determine all vertices a maximum of N edges
  /// away from each source vertex
  ///
  /// Templates:
  ///
  ///   Index     The type of the index into the face connectivity list
  ///
  /// Inputs:
  ///
  ///   N     The order of the neighborhood to calculate
  ///   F     #F by 3 face connectivity list
  ///
  /// Outputs:
  ///
  ///  NN     #V by 1 vector of neighborhoods for each vertex
  ///
  template <typename Index>
  NNI_INLINE void graphNeighborhood( const int N,
      const Eigen::Matrix<Index, Eigen::Dynamic, 3> &F,
      std::vector<std::vector<Index> > &NN );

}

#ifndef NNI_STATIC_LIBRARY
#  include "graphNeighborhood.cpp"
#endif

#endif
