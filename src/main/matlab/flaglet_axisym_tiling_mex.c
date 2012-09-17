// FLAGLET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include <flaglet.h>
#include <s2let.h>
#include "mex.h"

/**
 * MATLAB interface: flaglet_axisym_tiling_mex.
 * This function for internal use only.
 * Compute tiling in l-n harmonic space for axisymmetric wavelets.
 *
 * Usage: 
 *   [kappa kappa0] = flaglet_axisym_tiling_mex(B_l, B_n, L, N, J_min_l, J_min_n);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  int B_l, B_n, L, N, J_min_l, J_min_n;
  int iin = 0, iout = 0;

  // Check number of arguments
  if(nrhs!=6) {
    mexErrMsgIdAndTxt("flaglet_axisym_tiling_mex:InvalidInput:nrhs",
		      "Require six inputs.");
  }
  if(nlhs!=2) {
    mexErrMsgIdAndTxt("flaglet_axisym_tiling_mex:InvalidOutput:nlhs",
		      "Require two outputs.");
  }

  // Parse angular wavelet parameter B_l
  iin = 0;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_axisym_tiling_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B_l must be integer.");
  }
  B_l = (int)mxGetScalar(prhs[iin]);

  if (mxGetScalar(prhs[iin]) > (double)B_l || B_l <= 1)
    mexErrMsgIdAndTxt("flaglet_axisym_tiling_mex:InvalidInput:bandLimitNonInt",
          "Wavelet parameter B_l must be positive integer greater than 2");

  // Parse radial wavelet parameter B_n
  iin = 1;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_axisym_tiling_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B_n must be integer.");
  }
  B_n = (int)mxGetScalar(prhs[iin]);

  if (mxGetScalar(prhs[iin]) > (double)B_n || B_n <= 1)
    mexErrMsgIdAndTxt("flaglet_axisym_tiling_mex:InvalidInput:bandLimitNonInt",
          "Wavelet parameter B_n must be positive integer greater than 2");

  // Parse harmonic band-limit L
  iin = 2;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_axisym_tiling_mex:InvalidInput:LbandLimit",
		      "Harmonic band-limit L must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);


  if (mxGetScalar(prhs[iin]) > (double)L || L <= 0)
    mexErrMsgIdAndTxt("flaglet_axisym_tiling_mex:InvalidInput:bandLimitNonInt",
		      "Harmonic band-limit L must be positive integer.");

  if( B_l >= L ) {
    mexErrMsgIdAndTxt("flaglet_axisym_tiling_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B_l must be smaller than L!");
  }

    // Parse harmonic band-limit N
  iin = 3;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_axisym_tiling_mex:InvalidInput:LbandLimit",
          "Harmonic band-limit N must be integer.");
  }
  N = (int)mxGetScalar(prhs[iin]);


  if (mxGetScalar(prhs[iin]) > (double)N || N <= 0)
    mexErrMsgIdAndTxt("flaglet_axisym_tiling_mex:InvalidInput:bandLimitNonInt",
          "Harmonic band-limit N must be positive integer.");

  if( B_n >= N ) {
    mexErrMsgIdAndTxt("flaglet_axisym_tiling_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B_n must be smaller than N!");
  }
 
  // Parse angular first scale J_min_l
  iin = 4;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_axisym_tiling_mex:InvalidInput:J_min_l",
          "First scale J_min_l must be integer.");
  }
  J_min_l = (int)mxGetScalar(prhs[iin]);

  if (mxGetScalar(prhs[iin]) > (double)J_min_l || J_min_l < 0)
    mexErrMsgIdAndTxt("flaglet_axisym_tiling_mex:InvalidInput:J_min_l",
          "First scale J_min_l must be positive integer.");

  // Compute ultimate scale J_max
  int J_l = s2let_j_max(L, B_l);


  if( J_min_l > J_l+1 ) {
    mexErrMsgIdAndTxt("flaglet_axisym_tiling_mex:InvalidInput:J_min_l",
          "First scale J_min_l must be larger than that!");
  }

  // Parse angular first scale J_min_n
  iin = 5;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_axisym_tiling_mex:InvalidInput:J_min_n",
          "First scale J_min_n must be integer.");
  }
  J_min_n = (int)mxGetScalar(prhs[iin]);

  if (mxGetScalar(prhs[iin]) > (double)J_min_n || J_min_n < 0)
    mexErrMsgIdAndTxt("flaglet_axisym_tiling_mex:InvalidInput:J_min_n",
          "First scale J_min_n must be positive integer.");

  // Compute ultimate scale J_max
  int J_n = s2let_j_max(N, B_n);


  if( J_min_n > J_n+1 ) {
    mexErrMsgIdAndTxt("flaglet_axisym_tiling_mex:InvalidInput:J_min_n",
          "First scale J_min_n must be larger than that!");
  }

  // Allocate arrays
  double *kappa = (double*)calloc( (J_n+1) * (J_l+1) * L * N, sizeof(double));
  double *kappa0 = (double*)calloc( L * N, sizeof(double));

  // Run S2LET function
  flaglet_axisym_tiling(kappa, kappa0, B_l, B_n, L, N, J_min_l, J_min_n);

  // Output kappa and kappa0
  double *kappa_out, *kappa0_out;
  int n, l, jl, jn;

  iout = 0;
  plhs[iout] = mxCreateDoubleMatrix(1, (J_l+1)*(J_n+1)*L*N, mxREAL);
  kappa_out = mxGetPr(plhs[iout]);
  for (jn = J_min_n; jn <= J_n; jn++)
    for (jl = J_min_l; jl <= J_l; jl++)
      for (n = 0; n < N; n++)
        for (l = 0; l < L; l++){
          kappa_out[ jn*(J_l+1)*L*N  + jl*L*N + n*L + l ] = 
            creal(kappa[ jn*(J_l+1)*L*N  + jl*L*N + n*L + l ]);}

  iout = 1;
  plhs[iout] = mxCreateDoubleMatrix(N, L, mxREAL);
  kappa0_out = mxGetPr(plhs[iout]);
  for(n=0; n<N; n++)
    for(l=0; l<L; l++)
      kappa0_out[ l * N + n ] = creal(kappa0[ n * L + l ]); 

  free(kappa);
  free(kappa0);

}
