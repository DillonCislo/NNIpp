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

#ifndef _NNI_PARAM_H_
#define _NNI_PARAM_H_

#include "../General/nniInline.h"
#include <Eigen/Core>

namespace NNIpp {

  ///
  /// An enumeration of the different types of ghost point
  /// construction methods
  ///
  enum GHOST_POINT_CONSTRUCTION {

    // Custom user supplied ghost points
    NNI_GHOST_POINTS_CUSTOM = 1,

    // The dense circumcircle of the data bounding box
    NNI_GHOST_POINTS_CIRCLE = 2,

    // Rotated edges of the basic Delaunay triangulation
    NNI_GHOST_POINTS_EDGE = 3

  };

  ///
  /// An enumeration of the different gradient generation methods
  ///
  enum GRADIENT_GENERATION {

    // Direct gradient fitting to a high-order natural neighborhood
    NNI_GRADIENTS_DIRECT = 1,

    // Iterative fitting method to the first-order natural neighborhood
    NNI_GRADIENTS_ITER = 2

  };

  ///
  /// A class containing parameters and attributes used to streamline
  /// the construction of the 'NaturalNeighborInterpolant' object
  ///
  /// Templates:
  ///   Scalar    Input type of scattered data points, values, and gradients
  ///
  template <typename Scalar>
  class NNIParam {

    private:

      typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
      typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;

    public:

      // The number of scattered data points supplied
      int numPoints;

      // The method used for determining the position of the ghost points
      // used for natural neighbor extrapolation. The choice of method
      // will also dictate the type of gradient estimation used
      int ghostMethod;

      // A set of user defined custom ghost point coordinates
      Eigen::Matrix<Scalar, Eigen::Dynamic, 2> customGhostPoints;

      // The edge length increase factor for edge-based ghost point construction
      Scalar GPe;

      // The radius increase factor of the ghost point circle from the
      // circumcircle of the bounding box
      Scalar GPr;

      // The number of ghost points to create using the dense circle method
      int GPn;

      // The method used for gradient generation
      int gradType;

      // The relative importance of the two least-squares sub-problems
      // in the iterative form of Sibson's method for gradient generation
      Scalar iterAlpha;

      // If true, discrete gradient generation progress will be displayed
      bool dispGrad;

      // If true, interpolation progress will be displayed
      bool dispInterp;

      // #P by #V matrix. Holds the analytic x-gradient data of the input points
      Matrix DataGradX;

      // #P by #V matrix. Holds the analytic y-gradient data of the input points
      Matrix DataGradY;

      // Indicates whether analyic function gradients were supplied
      bool gradientSupplied;

      // #P by #V matrix. Holds the analytic xx-Hessian data of the input points
      Matrix DataHessXX;

      // #P by #V matrix. Holds the analytic xy-Hessian data of the input points
      Matrix DataHessXY;

      // #P by #V matrix. Holds the analytic yy-Hessian data of the input points
      Matrix DataHessYY;

      // Indicates whether analytic function Hessians were supplied
      bool hessianSupplied;

    public:

      ///
      /// Constructor for interpolation parameter class
      ///
      NNIParam(int nPts);

      ///
      /// Check the validity of the interpolation parameters
      ///
      NNI_INLINE void checkParam() const;

  };

}

#ifndef NNI_STATIC_LIBRARY
#  include "NNIParam.cpp"
#endif

#endif

