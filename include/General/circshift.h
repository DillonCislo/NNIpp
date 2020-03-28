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

#ifndef _CIRCSHIFT_H_
#define _CIRCSHIFT_H_

#include <cassert>
#include <cmath>
#include <iostream>

#include <Eigen/Core>

namespace NNIpp {

  ///
  /// Circularly shifts the values of a dense Eigen matrix/array.  Intended to mimic the
  /// functionality of MATLAB's `circshift' for 2D arrays.
  /// NOTE: If the row shift AND column shift are both non-zero then this algorithm
  /// shifts the rows FIRST and then the columns.
  ///
  /// Templates:
  ///   Derived   The derived type of the Eigen matrix/array (e.g. derived from MatrixXd)
  ///   Index     The derived type of the index variable (e.g. int )
  ///
  /// Inputs:
  ///   A         M by N Eigen matrix/array
  ///   a         An index variable that determines the row shift
  ///   b         An index variable that determines the column shift
  ///
  /// Outputs:
  ///   B         The shifted version of A
  ///
  template <typename Derived, typename Index>
  void circshift(
      Eigen::DenseBase<Derived> &B,
      const Eigen::DenseBase<Derived> &A,
      Index a, Index b = 0 ) {

    // Check that the input arguments are properly sized
    Index rows = A.rows();
    Index cols = A.cols();
    assert( ( rows == B.rows() ) && ( cols == B.cols() ) );

    // Find the sign of the shift variables
    int sgnA = ( Index(0) < a ) - ( a < Index(0) );
    int sgnB = ( Index(0) < b ) - ( b < Index(0) );

    // Bring the shift variables into the appropriate range
    Index aa = std::abs( a ) % rows;
    Index bb = std::abs( b ) % cols;

    // A temporary matrix/array object
    typename Eigen::DenseBase<Derived>::PlainObject rShift;
    if ( rShift.RowsAtCompileTime == -1 ) { rShift.resize( rows, Eigen::NoChange ); }
    if ( rShift.ColsAtCompileTime == -1 ) { rShift.resize( Eigen::NoChange, cols ); }


    // Perform the row shift
    switch( sgnA ) {

      case 1 :  // Move rows down

        rShift.topRows( aa ) = A.bottomRows( aa );
        rShift.bottomRows( rows - aa ) = A.topRows( rows - aa );
        break;

      case -1 : // Move rows up

        rShift.topRows( rows - aa ) = A.bottomRows( rows - aa );
        rShift.bottomRows( aa ) = A.topRows( aa );
        break;

      case 0 : // Do nothing

        rShift = A;
        break;

    }

    // Perform the column shift
    switch( sgnB ) {

      case 1 : // Move the columns right

        B.leftCols( bb ) = rShift.rightCols( bb );
        B.rightCols( cols - bb ) = rShift.leftCols( cols - bb );
        break;

      case -1 : // Move rows left

        B.leftCols( cols - bb ) = rShift.rightCols( cols - bb );
        B.rightCols( bb ) = rShift.leftCols( bb );
        break;

      case 0 : // Do nothing

        B = rShift;
        break;

    }

  };

} // namespace NNIpp

#endif // _CIRCSHIFT_H_
