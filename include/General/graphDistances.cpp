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

#include "graphDistances.h"

#include <cassert>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <limits>
#include <list>

#include <igl/adjacency_list.h>

template <typename DerivedF, typename DerivedD>
NNI_INLINE void NNIpp::graphDistances(
    const Eigen::MatrixBase<DerivedF> &F,
    Eigen::MatrixBase<DerivedD> &D ) {

  typedef typename Eigen::internal::traits<DerivedF>::Scalar Index;
  typedef typename Eigen::internal::traits<DerivedD>::Scalar Scalar;

  typedef Eigen::Array<bool, Eigen::Dynamic, 1> VectorXb;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;

  Scalar INF = std::numeric_limits<Scalar>::infinity();

  // The number of vertices
  Index numV = F.maxCoeff();
  numV += Index(1);

  // Construct the vertex adjacency list
  std::vector<std::vector<Index> > A(numV);
  igl::adjacency_list( F, A );

  D.resize( numV, numV );

  // Run a standard BFS for each vertex
  for( Index i = 0; i < numV; i++ ) {

    VectorXb visited = VectorXb::Constant( numV, 1, false );
    Vector d = Vector::Constant( numV, 1, INF );
    std::list<Index> queue;

    visited(i) = true;
    d(i) = Scalar(0);
    queue.push_back(i);

    while (!queue.empty()) {

      Index u = queue.front();
      queue.pop_front();

      for( Index j = 0; j < A[u].size(); j++ ) {

        if ( visited(A[u][j]) == false ) {

          visited(A[u][j]) = true;
          d(A[u][j]) = d(u) + Index(1);
          queue.push_back(A[u][j]);

        }

      }

    }

    D.col(i) = d;

  }

};

template <typename DerivedV, typename DerivedF, typename DerivedD>
NNI_INLINE void NNIpp::graphDistances(
    const Eigen::MatrixBase<DerivedV> &V,
    const Eigen::MatrixBase<DerivedF> &F,
    Eigen::MatrixBase<DerivedD> &D ) {

  std::runtime_error("This functionality has not been implemented yet!");

};

// TODO: Add explicit template instantiation
#ifdef NNI_STATIC_LIBRARY
#endif

