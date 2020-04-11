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

#include <Eigen/Core>

namespace NNIpp {

  ///
  /// Circularly shifts the values of a dense Eigen matrix/array.  Intended to mimic the
  /// functionality of MATLAB's `circshift' for 2D arrays.
  /// NOTE: If the row shift AND column shift are both non-zero then this algorithm
  /// shifts the rows FIRST and then the columns.
  ///
  /// Templates:
  ///
  ///   Derived   The derived type of the Eigen matrix/array (e.g. derived from MatrixXd)
  ///   Index     The derived type of the index variable (e.g. int )
  ///
  /// Inputs:
  ///
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
      Index a, Index b = 0 );

}

#endif
