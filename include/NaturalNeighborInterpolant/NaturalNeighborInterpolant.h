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

#ifndef _NATURAL_NEIGHBOR_INTERPOLANT_H_
#define _NATURAL_NEIGHBOR_INTERPOLANT_H_

#include <cassert>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <pair>

#include <Eigen/Core>

namespace NNIpp {

  ///
  /// An implementation of Sibson's C^1 continuous interpolant for scattered data
  /// using natural neighbor coordinates.  Supports C^1 continuous extrapolation
  /// using the assigned-value 'Ghost Point' method.  Output can be chosen to include
  /// both interpolated values and analytic derivatives of interpolants.  Supports
  /// multivariate interpolation.
  ///
  /// Templates:
  ///   Scalar    Input type of scattered data points, values, and gradients
  ///
  template <typename Scalar>
  class NaturalNeighborInterpolant {

    private:

      typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
      typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;

      // (#P+#GP) by #V matrix. Holds the function values of the scattered
      // data points and the ghost points
      Matrix m_Values;

      // (#P+#GP) by 2 matrix. Holds the scattered data and ghost points locations
      Eigen::Matrix<Scalar, Eigen::Dynamic, 2> m_Points;

      // #F by 3 matrix. Holds the face connectivity list of the extended
      // Delaunay triangulation
      Eigen::Matrix<int, Eigen::Dynamic, 3> m_Faces;

      // #E by 2 matrix. Holds the edge connectivity list of the extended
      // Delaunay triangulation
      Eigen::Matrix<int, Eigen::Dynamic, 2> m_Edges;

      // (#P+#GP) by #V matrix. Holds the x-gradient data of the input and ghost points
      Matrix m_DataGradX;

      // (#P+#GP) by matrix. Holds the y-gradient data of the input and ghost points
      Matrix m_DataGradY;

      // #F by 1 vector. Holds the 'Delta' parameter of each face of the extended
      // Delaunay triangulation.  Equal to twice the signed area of each face
      Vector m_fDelta;

      // #GP by 2 matrix.  Holds the coordinates of the convex hull of the extended
      // Delaunay triangulation
      Eigen::Matrix<Scalar, Eigen::Dynamic, 2> m_ConvexHull;

      // Display progress towards estimating discrete gradients
      bool m_dispGrad;

      // Display progress towards interpolating values
      bool m_dispInterp;


    public:

      ///
      /// Constructor with user supplied analytic gradient data
      ///
      /// Inputs:
      ///
      ///   Xp            #P by 1 list of data point x-coordinates
      ///   Yp            #P by 1 list of data point y-coordinates
      ///   Vp            #P by #V matrix of data point function values
      ///   DVx           #P by #V matrix of analytic x-gradient values at data points
      ///   DVy           #P by #V matrix of analytic y-gradient values at data points
      ///   GPn           The number of ghost points to create
      ///   GPr           The radius increase factor of the ghost point circle
      ///                 from the circumcircle of the data point bounding box
      ///   dispInterp    Boolean display interpolation progress
      ///
      NaturalNeighborInterpolant( const Vector &Xp, const Vector &Yp,
          const Matrix &Vp, const Matrix &DVx, const Matrix &DVy,
          int GPn, Scalar GPr, bool dispInterp );

      ///
      /// Constructor with discrete gradient generation
      ///
      /// Inputs:
      ///
      ///   Xp            #P by 1 list of data point x-coordinates
      ///   Yp            #P by 1 list of data point y-coordinates
      ///   Vp            #P by #V matrix of data point function values
      ///   DVx           #P by #V matrix of analytic x-gradient values at data points
      ///   DVy           #P by #V matrix of analytic y-gradient values at data points
      ///   GPn           The number of ghost points to create
      ///   GPr           The radius increase factor of the ghost point circle
      ///                 from the circumcircle of the data point bounding box
      ///   dispGrad      Boolean display gradient generation progress
      ///   dispInterp    Boolean display interpolation progress
      ///
      NaturalNeighborInterpolant( const Vector &Xp, const Vector &Yp,
          const Matrix &Vp, int GPn, Scalar GPr,
          bool dispGrad, bool dispInterp );

      /// Interpolate function at a set of query points using Sibson's
      /// C^1 continuous method for scattered data interpolation
      ///
      /// Inputs:
      ///
      ///   Xq      #Q by 1 list of query point x-coordinates
      ///   Yq      #Q by 1 list of query point y-coordinates
      ///
      /// Outputs:
      ///
      ///   Fq      #Q by #V matrix of interpolated function values
      ///
      void operator()( const Vector &Xq, const Vector &Yq, Matrix &Fq );

      /// Interpolate function at a set of query points using Sibson's
      /// C^1 continuous method for scattered data interpolation.
      /// Also evaluate analytic function gradients
      ///
      /// Inputs:
      ///
      ///   Xq      #Q by 1 list of query point x-coordinates
      ///   Yq      #Q by 1 list of query point y-coordinates
      ///
      /// Outputs:
      ///
      ///   Fq      #Q by #V matrix of interpolated function values
      ///   DFx     #Q by #V matrix of x-derivatives of funciton values
      ///   DFy     #Q by #V matrix of y-derivatives of function values
      ///
      void operator()( const Vector &Xq, const Vector &Yq,
          Matrix &Fq, Matrix &DFx, Matrix &DFy );

      ///
      /// Calculate the natural neighbor coordinates of a set of query points
      ///
      /// Inputs:
      ///
      ///   Xq      #Q by 1 list of query point x-coordinates
      ///   Yq      #Q by 1 list of query point y-coordinates
      ///
      /// Outputs:
      ///
      ///   u       #Q by 1 vector. Each entry is a list of the (normalized)
      ///           natural neighbor coordinates of the query point
      ///   uIDx    #Q by 1 vector. Each entry is a list of the vertex IDs
      ///           of the corresponding natural neighbor point in the
      ///           extended Delaunay triangulation
      ///
      void naturalNeighborCoordinates( const Vector &Xq, const Vector &Yq,
          std::vector<Vector> &u, std::vector<Eigen::VectorXi> &uIDx );

      ///
      /// Calculate the natural neighbor coordinates of a set of query points
      /// Also extract the sufficient geometric data to calculate the derivative
      /// of the natural neighbor coordinates with respect to the query point
      /// coordinates
      ///
      /// Inputs:
      ///
      ///   Xq      #Q by 1 list of query point x-coordinates
      ///   Yq      #Q by 1 list of query point y-coordinates
      ///
      /// Outputs:
      ///
      ///   u       #Q by 1 vector. Each entry is a list of the (normalized)
      ///           natural neighbor coordinates of the query point
      ///   uIDx    #Q by 1 vector. Each entry is a list of the vertex IDs
      ///           of the corresponding natural neighbor point in the
      ///           extended Delaunay triangulation
      ///   uVC     #Q by 1 vector. Each entry is another a matrix whose row
      ///           entries define the counter-clockwise ordered vertices of
      ///           the virtual Voronoi cell of the corresponding query point
      ///   uA      #Q by 1 list of the areas of the virtual Voronoi cells
      ///           of each query point
      ///
      void naturalNeighborCoordinates( const Vector &Xq, const Vector &Yq,
          std::vector<Vector> &u, std::vector<Eigen::VectorXi> &uIDx,
          std::vector<Eigen::Matrix<Scalar, Eigen::Dynamic, 2> > &uVC,
          Vector &uA );

    private:

      ///
      /// Generate discrete gradients. Gradients at a data point are defined
      /// to be the first-order coefficients of a third-order Taylor polynomial
      /// fit to the third-degree natural neighborhood of that same data point.
      ///
      /// Inputs:
      ///
      ///   numPoints   The number of data points supplied (excludes ghost points)
      ///
      void generateGradients( int numPoints );

  };

} // namespace NNIpp

#endif // _NATURAL_NEIGHBOR_INTERPOLANT_H_
