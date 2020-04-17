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

#include "graphNeighborhood.h"

#include <cassert>
#include <stdexcept>
#include <iostream>
#include <list>

#include <igl/adjacency_list.h>

template <typename Index>
NNI_INLINE void NNIpp::graphNeighborhood( const int N,
    const Eigen::Matrix<Index, Eigen::Dynamic, 3> &F,
    std::vector<std::vector<Index> > &NN ) {

  assert( (N > Index(0)) && "Neighborhood order must be positive" );

  typedef Eigen::Array<bool, Eigen::Dynamic, 1> VectorXb;
  typedef Eigen::Matrix<Index, Eigen::Dynamic, 1> Vector;

  // The number of vertices
  Index numV = F.maxCoeff() + Index(1);

  // Construct the basic vertex adjacency list
  std::vector<std::vector<Index> > A(numV);
  igl::adjacency_list(F, A);

  if ( N == 1 ) {

    NN = A;
    return;

  }

  NN.resize(numV);

  // Run a standard BFS for each vertex terminating at order N
  for( Index i = 0; i < numV; i++ ) {

    std::vector<Index> curNN;

    VectorXb visited = VectorXb::Constant( numV, 1, false );
    Vector d = Vector::Constant( numV, 1, Index(-1) );
    std::list<Index> queue;

    visited(i) = true;
    d(i) = Index(0);
    queue.push_back(i);

    while(!queue.empty()) {

      Index u = queue.front();

      if ( d(u) > N ) {
        break;
      } else {
        curNN.push_back(u);
      }

      queue.pop_front();

      for( Index j = 0; j < A[u].size(); j++ ) {

        if ( visited(A[u][j]) == false ) {

          visited(A[u][j]) = true;
          d(A[u][j]) = d(u) + Index(1);
          queue.push_back(A[u][j]);

        }

      }

    }

    // Remove self-reference in neighborhood
    curNN.erase( curNN.begin() );
    NN[i] = curNN;

  }

};

// TODO: Add explicit template instantiation
#ifdef NNI_STATIC_LIBRARY
#endif
