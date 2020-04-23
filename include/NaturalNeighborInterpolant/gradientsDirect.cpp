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

#include "gradientsDirect.h"

#include <cassert>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "../General/graphNeighborhood.h"

template <typename Scalar>
NNI_INLINE void NNIpp::gradientsDirect(
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 2> &V,
    const Eigen::Matrix<int, Eigen::Dynamic, 3> &F,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &Vals,
    const NNIParam<Scalar> &param,
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &Dx,
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &Dy,
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &Dxx,
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &Dxy,
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &Dyy ) {

  typedef Eigen::Matrix<Scalar, 1, 2> Vector2d;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
  typedef Eigen::Array<Scalar, 1, Eigen::Dynamic> ArrayRow;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;

  // The number of points in the extended triangulation
  int numAllPoints = V.rows();
  
  // The number of input points
  int numPoints = param.numPoints;

  // The number of values per point
  int numVals = param.numVals;

  // Calculate the 3rd-order natural neighboorhood of each input point
  std::vector<std::vector<int> > NN( numAllPoints );
  NNIpp::graphNeighborhood( 3, F, NN );

  // The inverse weights of the Taylor polynomial
  ArrayRow a(9);
  a << Scalar(1.0), Scalar(1.0), Scalar(2.0), Scalar(1.0),
    Scalar(2.0), Scalar(3.0), Scalar(2.0), Scalar(2.0), Scalar(3.0);
  a = a.inverse();

  // Iterate over all input points
  for( int i = 0; i < numPoints; i++ ) {

    // The current point
    Vector2d curV = V.row(i);

    // The natural neighborhood of the current point
    std::vector<int> curNN = NN[i];
    
    // Remove ghost points from consideration
    std::vector<int>::iterator it;
    it = std::remove_if( curNN.begin(), curNN.end(),
        std::bind2nd( std::greater<int>(), (numPoints-1) ) );
    curNN.erase( it, curNN.end() );

    // The number of natural neighbors for the current point
    int numNN = curNN.size();

    // Construct the least squares LHS matrix
    Matrix A(numNN, 9);
    Vector invL(numNN);
    for( int j = 0; j < numNN; j++ ) {

      // Components of the separation vectors
      Scalar X = V(curNN[j], 0) - curV(0);
      Scalar Y = V(curNN[j], 1) - curV(1);

      // Inverse separation length
      invL(j) = Scalar(1.0) / std::sqrt( X*X + Y*Y );

      ArrayRow Arow(9);
      Arow << X, Y, X*X, X*Y, Y*Y, X*X*X, X*X*Y, X*Y*Y, Y*Y*Y;
      Arow = invL(j) * a * Arow;

      A.row(j) = Arow.matrix();

    }

    // Solve the linear system for each value component
    for( int j = 0; j < numVals; j++ ) {

      // The least-squares RHS vector
      Vector B(numNN);
      for( int k = 0; k < numNN; k++ ) {
        B(k) = invL(k) * ( Vals(curNN[k], j) - Vals(i,j) );
      }

      // Solve the linear least-squares problem
      Eigen::Matrix<Scalar, 9, 1> g;

      // LEAST-SQUARES METHOD 1: SVD DECOMPOSITION
      g = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(B);

      // LEAST-SQUARES METHOD 2: QR NO PIVOTING
      // g = A.householderQr().solve(B);

      // LEAST-SQUARES METHOD 3: QR COLUMN PIVOTING
      // g = A.colPivHouseholderQr().solve(B);

      // LEAST-SQUARES METHOD 4: QR FULL PIVOTING
      // g = A.fullPivHouseholderQr().solve(B);

      // LEAST-SQUARES METHOD 5: NORMAL EQUATIONS
      // g = ( A.transpose() * A ).ldlt().solve( A.transpose() * B );
      
      Dx(i,j) = g(0);
      Dy(i,j) = g(1);

      Dxx(i,j) = g(2);
      Dxy(i,j) = g(3);
      Dyy(i,j) = g(4);

    }

  }

};
      
// TODO: Add explicit template instantiation
#ifdef NNI_STATIC_LIBRARY
#endif
