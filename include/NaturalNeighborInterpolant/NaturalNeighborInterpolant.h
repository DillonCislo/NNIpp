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

#include "../General/nniInline.h"

#include <vector>
#include <Eigen/Core>
#include "NNIParam.h"
#include <igl/AABB.h>

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
      typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> ArrayVec;
      typedef Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> Array;

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

      // #E by 2 list of facets, such that m_OrientedEdges.row(f+#F*c)
      // is the edge opposite F(f,c)
      Eigen::Matrix<int, Eigen::Dynamic, 2> m_OrientedEdges;

      // #F by 3 face adjacency list. m_FaceNeighbors(i,j) is the ID of
      // the face adjacent to the jth edge of face i, i.e. the face
      // opposite vertex j
      Eigen::Matrix<int, Eigen::Dynamic, 3> m_FaceNeighbors;

      // An AABB tree. Used to determine which mesh face contains a query point
      igl::AABB<Matrix, 2> m_Tree;

      // #F by 3 matrix. Holds the x-coordinates of each vertex in each face
      // of the extended triangulation
      Eigen::Matrix<Scalar, Eigen::Dynamic, 3> m_FaceVertexX;

      // #F by 3 matrix. Holds the y-coordinates of each vertex in each face
      // of the extended triangulation
      Eigen::Matrix<Scalar, Eigen::Dynamic, 3> m_FaceVertexY;

      // #F by 2 matrix. Holds the (x,y)-coordinates of the circumcenter
      // of each face of the extended triangulation
      Eigen::Matrix<Scalar, Eigen::Dynamic, 2> m_FCC;

      // #F by 1 vector. Holds the 'Delta' parameter of each face of the extended
      // Delaunay triangulation.  Equal to twice the signed area of each face
      ArrayVec m_fDelta;

      // #F by 4 matrix. Holds precomputed data needed to calculate the 'Gamma'
      // parameter for each face of the extended triangulation given a query point
      Eigen::Array<Scalar, Eigen::Dynamic, 4> m_fGamma;

      // #CP by 2 matrix.  Holds the coordinates of the convex hull of the extended
      // Delaunay triangulation
      Eigen::Matrix<Scalar, Eigen::Dynamic, 2> m_ConvexHull;

      // (#P+#GP) by #V matrix. Holds the x-gradient data of the input and ghost points
      Matrix m_DataGradX;

      // (#P+#GP) by #V matrix. Holds the y-gradient data of the input and ghost points
      Matrix m_DataGradY;

      // (#P+#GP) by #V matrix. Holds the xx-Hessian data of the input and ghost points
      Matrix m_DataHessXX;

      // (#P+#GP) by #V matrix. Holds the xy-Hessian data of the input and ghost points
      Matrix m_DataHessXY;

      // (#P+#GP) by #V matrix. Holds the yy-Hessian data of the input and ghost points
      Matrix m_DataHessYY;

      // Display progress towards estimating discrete gradients
      bool m_dispGrad;

      // Display progress towards interpolating values
      bool m_dispInterp;

  public:

      ///
      /// Null constructor
      ///
      NaturalNeighborInterpolant() {};

      ///
      /// Default constructor
      ///
      /// Inputs:
      ///
      ///   Xp            #P by 1 list of data point x-coordinates
      ///   Yp            #P by 1 list of data point y-coordinates
      ///   Vp            #P by #V matrix of data point function values
      ///   param         An 'NNIParam' class containing the rest of the
      ///                 parameters needed to construct the interpolant
      ///
      NaturalNeighborInterpolant( const Vector &Xp, const Vector &Yp,
          const Matrix &Vp, const NNIParam<Scalar> &param );

      ///
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
      void operator()( const Vector &Xq, const Vector &Yq, Matrix &Fq ) const;

      ///
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
          Matrix &Fq, Matrix &DFx, Matrix &DFy ) const;

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
          std::vector<Vector> &u, std::vector<Eigen::VectorXi> &uIDx ) const;

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
      ///           shifted so that the query point is at the origin
      ///   uA      #Q by 1 list of the areas of the virtual Voronoi cells
      ///           of each query point
      ///
      void naturalNeighborCoordinates( const Vector &Xq, const Vector &Yq,
          std::vector<Vector> &u, std::vector<Eigen::VectorXi> &uIDx,
          std::vector<Eigen::Matrix<Scalar, Eigen::Dynamic, 2> > &uVC,
          Vector &uA ) const;

    private:

      ///
      /// Generate ghost points
      ///
      /// Inputs:
      ///
      ///   Xp            #P by 1 list of data point x-coordinates
      ///   Yp            #P by 1 list of data point y-coordinates
      ///   param         An 'NNIParam' class containing the rest of the
      ///                 parameters needed to construct the interpolant
      ///
      /// Outputs:
      ///
      ///   GPx           #GP by 1 list of ghost point x-coordinates
      ///   GPy           #GP by 1 list of ghost point y-coordinates
      ///
      NNI_INLINE void generateGhostPoints(
          const Vector &Xp, const Vector &Yp,
          const NNIParam<Scalar> &param,
          Vector &GPx, Vector &GPy );

      ///
      /// Generate discrete gradients
      ///
      /// Inputs:
      ///
      ///   param         An 'NNIParam' class containing the rest of the
      ///                 parameters needed to construct the interpolant
      ///
      NNI_INLINE void generateGradients( const NNIParam<Scalar> &param );

      ///
      /// Generate values/derivatives for ghost points
      ///
      /// Inputs:
      ///
      ///   param         An 'NNIParam' class containing the rest of the
      ///                 parameters needed to construct the interpolant
      ///
      NNI_INLINE void ghostPointValueHandling( const NNIParam<Scalar> &param );

      ///
      /// Precompute data for calculating the 'Gamma'
      /// parameter on each face of the extended triangulation
      ///
      NNI_INLINE void precomputeGamma();

      ///
      /// Calculate the 'Gamma' parameter from (Hiyoshi, 2008). Used to
      /// calculate the natural neighbor coordinates of a query point
      /// Gamma(v1, v2, v3, v4) does not have a simple geometric interpretation
      ///
      /// Inputs:
      ///
      ///   X   X-coordinate of query point
      ///   Y   Y-coordinate of query point
      ///
      /// Outputs:
      ///
      ///   G   #F by 1 list of 'Gamma' values
      ///
      NNI_INLINE void gamma( const Scalar X, const Scalar Y, ArrayVec &G ) const;

      ///
      /// Determine if a query point lies within the circumcircle of a face
      /// of the extended triangulation
      ///
      /// Inputs:
      ///
      ///   Xq    X-coordinate of query point
      ///   Yq    Y-coordinate of query point
      ///   FID   The ID of the face being checked
      ///
      /// Outputs:
      ///
      ///   in    True if (Xq,Yq) lies within the circumcircle
      ///
      NNI_INLINE bool inCircle( const Scalar Xq, const Scalar Yq, const int FID ) const;

      ///
      /// Determine whether each point in a list of query points lies within
      /// the faces of the extended triangulation
      ///
      /// Inputs:
      ///
      ///   X   #Q by 1 list of query point x-coordinates
      ///   Y   #Q by 1 list of query point y-coordinates
      ///
      /// Outputs:
      ///
      ///   I   #Q by 1 list of indices into the face connectivity list
      ///       of the first containing element (-1 means no containing
      ///       element)
      ///
      NNI_INLINE void inElement(
          const Vector &Xq, const Vector &Yq, Eigen::VectorXi &I ) const;

  };

}

#ifndef NNI_STATIC_LIBRARY
#  include "NaturalNeighborInterpolant.cpp"
#endif

#endif
