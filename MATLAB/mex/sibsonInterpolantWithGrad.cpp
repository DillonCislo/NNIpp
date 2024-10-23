/* =============================================================================================
 *
 *  sibsonInterpolantWithGrad.cpp
 *  
 *
 *
 *  by Dillon Cislo
 *  04/25/2020
 *
 *  This is a MEX-file for MATLAB
 *  
 * ============================================================================================*/

#include "mex.h" // for MATLAB

#include <Eigen/Core>
#include <chrono>
#include <iostream>

#include "../../include/NaturalNeighborInterpolant/NaturalNeighborInterpolant.h"
#include "../../include/NaturalNeighborInterpolant/NNIParam.h"

// Main function
void mexFunction( int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[] ) {

  // Check for proper number of arguments
  if ( nrhs != 6 ) {
    mexErrMsgIdAndTxt( "MATLAB:sibsoninterpolantwithgrad:nargin",
        "SIBSONINTERPOLANTWITHGRAD requires 6 input arguments" );
  } else if ( nlhs != 3 ) {
    mexErrMsgIdAndTxt("MATLAB:sibsoninterpolantwithgrad:nargout",
        "SIBSONINTERPOLANTWITHGRAD requries 3 output arguments" );
  }

  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 2> PointVector;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;

  // X-coordinates of input data points
  double *xIn = mxGetPr( prhs[0] );
  int numPoints  = (int) mxGetM( prhs[0] );

  // Y-coordinates of input data points
  double *yIn = mxGetPr( prhs[1] );

  // The values to be interpolated
  double *vIn = mxGetPr( prhs[2] );
  int numVals = (int) mxGetN( prhs[2] );

  // X-coordinates of query points
  double *xQin = mxGetPr( prhs[4] );
  int numQueries = (int) mxGetM( prhs[4] );

  // Y-coordinates of query points
  double *yQin = mxGetPr( prhs[5] );

  // Map input arrays to Eigen-style matrices
  Vector Xp = Eigen::Map<Vector>(xIn, numPoints, 1);
  Vector Yp = Eigen::Map<Vector>(yIn, numPoints, 1);
  Matrix Vp = Eigen::Map<Matrix>(vIn, numPoints, numVals);
  Vector Xq = Eigen::Map<Vector>(xQin, numQueries, 1);
  Vector Yq = Eigen::Map<Vector>(yQin, numQueries, 1);

  // Parse input options
  NNIpp::NNIParam<double> param( numPoints, numVals );

  int idx, tmp;

  // The method used for constructing ghost points
  if ( (idx = mxGetFieldNumber( prhs[3], "ghostMethod" )) != -1 ) {

    tmp = (int) *mxGetPr(mxGetFieldByNumber( prhs[3], 0, idx ));

    if (tmp == 1) {
      param.ghostMethod = NNIpp::NNI_GHOST_POINTS_CUSTOM;
    } else if (tmp == 2) {
      param.ghostMethod = NNIpp::NNI_GHOST_POINTS_CIRCLE;
    } else if (tmp == 3) {
      param.ghostMethod = NNIpp::NNI_GHOST_POINTS_EDGE;
    } else {
      mexErrMsgTxt("Invalid ghost point construction method supplied!");
    }

  }

  // Custom ghost point coordinates
  if ( (idx = mxGetFieldNumber( prhs[3], "customGhostPoints" )) != -1 ) {

    int GPn = (int) mxGetM(mxGetFieldByNumber( prhs[3], 0, idx ));

    param.customGhostPoints =
      Eigen::Map<PointVector>( mxGetPr(mxGetFieldByNumber( prhs[3], 0, idx )), GPn, 1 );

    // Override other directives for ghost point construction
    param.ghostMethod = NNIpp::NNI_GHOST_POINTS_CUSTOM;

  }

  // Edge length increase factor used for edge-based ghost point construction
  if ( (idx = mxGetFieldNumber( prhs[3], "GPe" )) != -1 ) {
    param.GPe = *mxGetPr(mxGetFieldByNumber( prhs[3], 0, idx ));
  }

  // Radius increase factor of ghost point circle
  if ( (idx = mxGetFieldNumber( prhs[3], "GPr" )) != -1) {
    param.GPr = *mxGetPr(mxGetFieldByNumber( prhs[3], 0, idx ));
  }

  // The number of ghost points to create using the dense circle method
  if ( (idx = mxGetFieldNumber( prhs[3], "GPn" )) != -1) {
    param.GPn = (int) *mxGetPr(mxGetFieldByNumber( prhs[3], 0, idx ));
  }

  // The method used for gradient generation
  if ( (idx = mxGetFieldNumber( prhs[3], "gradType" )) != -1 ) {
    
    tmp = (int) *mxGetPr(mxGetFieldByNumber( prhs[3], 0, idx ));

    if (tmp == 1) {
      param.gradType = NNIpp::NNI_GRADIENTS_DIRECT;
    } else if (tmp == 2) {
      param.gradType = NNIpp::NNI_GRADIENTS_ITER;
    } else {
      mexErrMsgTxt("Invalid gradient generation method supplied!");
    }
  
  }

  // Parameter for Sibson's method for gradient generation
  if ( (idx = mxGetFieldNumber( prhs[3], "iterAlpha" )) != -1 ) {
    param.iterAlpha = *mxGetPr(mxGetFieldByNumber( prhs[3], 0, idx ));
  }

  // Display option for gradient generation progress
  if ( (idx = mxGetFieldNumber( prhs[3], "dispGrad" )) != -1 ) {
    param.dispGrad = *mxGetLogicals(mxGetFieldByNumber( prhs[3], 0, idx));
  }

  // Display option for interpolation progress
  if ( (idx= mxGetFieldNumber( prhs[3], "dispInterp" )) != -1 ) {
    param.dispInterp = *mxGetLogicals(mxGetFieldByNumber( prhs[3], 0, idx));
  }

  int idX = mxGetFieldNumber( prhs[3], "DataGradX" );
  bool hasDX = idX != -1;

  int idY = mxGetFieldNumber( prhs[3], "DataGradY" );
  bool hasDY = idY != -1;

  if ( hasDX && hasDY ) {

    param.DataGradX =
      Eigen::Map<Matrix>(mxGetPr(mxGetFieldByNumber(prhs[3], 0, idX)), numPoints, numVals);

    param.DataGradY =
      Eigen::Map<Matrix>(mxGetPr(mxGetFieldByNumber(prhs[3], 0, idY)), numPoints, numVals);

    param.gradientSupplied = true;

  } else if ( hasDX || hasDY ) {

    mexErrMsgTxt("Incomplete gradient data supplied");
  
  }

  int idXX = mxGetFieldNumber( prhs[3], "DataHessXX" );
  bool hasDXX = idXX != -1;

  int idXY = mxGetFieldNumber( prhs[3], "DataHessXY" );
  bool hasDXY = idXY != -1;

  int idYY = mxGetFieldNumber( prhs[3], "DataHessYY" );
  bool hasDYY = idYY != -1;

  if ( hasDXX && hasDXY && hasDYY ) {

    param.DataHessXX =
      Eigen::Map<Matrix>(mxGetPr(mxGetFieldByNumber(prhs[3], 0, idXX)), numPoints, numVals);

    param.DataHessXY =
      Eigen::Map<Matrix>(mxGetPr(mxGetFieldByNumber(prhs[3], 0, idXY)), numPoints, numVals);

    param.DataHessYY =
      Eigen::Map<Matrix>(mxGetPr(mxGetFieldByNumber(prhs[3], 0, idYY)), numPoints, numVals);

    param.hessianSupplied = true;

  } else if ( hasDXX || hasDYY || hasDYY ) {

    mexErrMsgTxt("Incomplete Hessian data supplied!");

  }

  // Generate an interpolant
  NNIpp::NaturalNeighborInterpolant<double> NNI( Xp, Yp, Vp, param );

  // Evaluate interpolant at query points
  Matrix Fq( numQueries, numVals );
  Matrix DFx( numQueries, numVals );
  Matrix DFy( numQueries, numVals );

  // auto start = std::chrono::high_resolution_clock::now();
  
  try {

    NNI( Xq, Yq, Fq, DFx, DFy );

  } catch (std::runtime_error& e) {

    // std::cerr << e.what() << std::endl;
    // return -1
    
    mexErrMsgTxt(e.what());

  }

  // auto stop = std::chrono::high_resolution_clock::now();
  // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);

  // std::cout << (1.0e-6 * duration.count()) << " seconds." << std::endl;

  // Output interpolated values
  plhs[0] = mxCreateDoubleMatrix( numQueries, numVals, mxREAL );
  Eigen::Map<Eigen::MatrixXd>(mxGetPr(plhs[0]), numQueries, numVals) = Fq;

  // Output x-derivatives
  plhs[1] = mxCreateDoubleMatrix( numQueries, numVals, mxREAL );
  Eigen::Map<Eigen::MatrixXd>(mxGetPr(plhs[1]), numQueries, numVals) = DFx;

  // Output y-derivatives
  plhs[2] = mxCreateDoubleMatrix( numQueries, numVals, mxREAL );
  Eigen::Map<Eigen::MatrixXd>(mxGetPr(plhs[2]), numQueries, numVals) = DFy;

  return;

};
