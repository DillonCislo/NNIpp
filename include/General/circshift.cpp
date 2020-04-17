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

#include "circshift.h"

#include <cassert>
#include <cmath>
#include <iostream>

template <typename Derived, typename Index>
NNI_INLINE void NNIpp::circshift(
    Eigen::DenseBase<Derived> &B,
    const Eigen::DenseBase<Derived> &A,
    Index a, Index b ) {

  Index rows = A.rows();
  Index cols = A.cols();
  B = Eigen::DenseBase<Derived>::Zero( rows, cols );

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

// TODO: Add explicit template instantiaion
#ifdef NNI_STATIC_LIBRARY
#endif
