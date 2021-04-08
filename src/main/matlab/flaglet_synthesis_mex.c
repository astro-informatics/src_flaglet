// FLAGLET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include <flaglet.h>
#include "mex.h"
#include <math.h>
#include <s2let.h>
#include <complex.h>

/**
 * MATLAB interface: flaglet_axisym_synthesis.
 * This function for internal use only.
 * Compute axisymmetric wavelet transform (synthesis)
 * with output in pixel space.
 *
 * Usage: 
 *   f = ...
 *        flaglet_axisym_synthesis_mex(f_wav, f_scal, B_l, B_p, L, P, J_min_l, J_min_p, N, Downsample);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  int n, p, t, i, j, f_m, f_n;
  int B_l, B_p, L, P, N, J_min_l, J_min_p;
  flaglet_parameters_t parameters = {};
  double *f_wav_real, *f_scal_real, *f_real, *f_wav_imag, *f_scal_imag, *f_imag;
  complex double *f_wav, *f_scal, *f;
  int iin = 0, iout = 0;

  // Check number of arguments
  if(nrhs!=11) {
    mexErrMsgIdAndTxt("flaglet_axisym_synthesis_mex:InvalidInput:nrhs",
          "Require eleven inputs.");
  }
  if(nlhs!=1) {
    mexErrMsgIdAndTxt("flaglet_axisym_synthesis_mex:InvalidOutput:nlhs",
          "Require two outputs.");
  }

  // Parse multiresolution flag
  iin = 10;
  if( !mxIsLogicalScalar(prhs[iin]) )
    mexErrMsgIdAndTxt("flaglet_axisym_synthesis_mex:InvalidInput:downsample",
          "Multiresolution flag must be logical.");
  parameters.upsample = !mxIsLogicalScalarTrue(prhs[iin]);

  // Parse input wavelets f_wav
  iin = 0;
  if( !mxIsDouble(prhs[iin]) ) {
    mexErrMsgIdAndTxt("flag_analysis_mex:InvalidInput:f",
          "Function values must be doubles.");
  }
  f_m = mxGetM(prhs[iin]);
  f_n = mxGetN(prhs[iin]);
  int f_is_complex = mxIsComplex(prhs[iin]);
  f_wav_real = mxGetPr(prhs[iin]);
  f_wav_imag = mxGetPi(prhs[iin]);
  f_wav = (complex double*)malloc( f_m*f_n * sizeof(complex double));
  for(j=0; j<f_m*f_n; j++)
      f_wav[ j ] = f_wav_real[ j ] + I * (f_is_complex ? f_wav_imag[ j ] : 0.0);

  // Parse input scaling function f_scal
  iin = 1;
  f_m = mxGetM(prhs[iin]);
  f_n = mxGetN(prhs[iin]);
  f_is_complex = mxIsComplex(prhs[iin]);
  f_scal_real = mxGetPr(prhs[iin]);
  f_scal_imag = mxGetPi(prhs[iin]);
  f_scal = (complex double*)malloc( f_m*f_n * sizeof(complex double));
  for(t=0; t<f_m; t++)
    for(p=0; p<f_n; p++)
      f_scal[ t*f_n + p ] = f_scal_real[ p*f_m + t ] + I * (f_is_complex ? f_scal_imag[ p*f_m + t ] : 0.0);

  // Parse angular wavelet parameter B_l
  iin = 2;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_axisym_synthesis_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B_l must be integer.");
  }
  B_l = (int)mxGetScalar(prhs[iin]);
  parameters.B_l = B_l;

  if (mxGetScalar(prhs[iin]) > (double)B_l || B_l <= 1)
    mexErrMsgIdAndTxt("flaglet_axisym_synthesis_mex:InvalidInput:bandLimitNonInt",
          "Wavelet parameter B_l must be positive integer greater than 2");

  // Parse radial wavelet parameter B_p
  iin = 3;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_axisym_synthesis_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B_p must be integer.");
  }
  B_p = (int)mxGetScalar(prhs[iin]);
  parameters.B_p = B_p;

  if (mxGetScalar(prhs[iin]) > (double)B_p || B_p <= 1)
    mexErrMsgIdAndTxt("flaglet_axisym_synthesis_mex:InvalidInput:bandLimitNonInt",
          "Wavelet parameter B_p must be positive integer greater than 2");

  // Parse harmonic band-limit L
  iin = 4;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_axisym_synthesis_mex:InvalidInput:LbandLimit",
          "Harmonic band-limit L must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);
  parameters.L = L;


  if (mxGetScalar(prhs[iin]) > (double)L || L <= 0)
    mexErrMsgIdAndTxt("flaglet_axisym_synthesis_mex:InvalidInput:bandLimitNonInt",
          "Harmonic band-limit L must be positive integer.");

  if( B_l >= L ) {
    mexErrMsgIdAndTxt("flaglet_axisym_synthesis_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B_l must be smaller than L!");
  }

    // Parse harmonic band-limit P
  iin = 5;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_axisym_synthesis_mex:InvalidInput:LbandLimit",
          "Harmonic band-limit P must be integer.");
  }
  P = (int)mxGetScalar(prhs[iin]);
  parameters.P = P;


  if (mxGetScalar(prhs[iin]) > (double)P|| P <= 0)
    mexErrMsgIdAndTxt("flaglet_axisym_synthesis_mex:InvalidInput:bandLimitNonInt",
          "Harmonic band-limit P must be positive integer.");

  if( B_p >= P ) {
    mexErrMsgIdAndTxt("flaglet_axisym_synthesis_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B_p must be smaller than P!");
  }
 
  // Parse angular first scale J_min_l
  iin = 6;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_axisym_synthesis_mex:InvalidInput:J_min_l",
          "First scale J_min_l must be integer.");
  }
  J_min_l = (int)mxGetScalar(prhs[iin]);
  parameters.J_min_l = J_min_l;

  if (mxGetScalar(prhs[iin]) > (double)J_min_l || J_min_l < 0)
    mexErrMsgIdAndTxt("flaglet_axisym_synthesis_mex:InvalidInput:J_min_l",
          "First scale J_min_l must be positive integer.");

  // Compute ultimate scale J_max
  int J_l = ceil(log(L) / log(B_l));


  if( J_min_l > J_l+1 ) {
    mexErrMsgIdAndTxt("flaglet_axisym_synthesis_mex:InvalidInput:J_min_l",
          "First scale J_min_l must be larger than that!");
  }

  // Parse angular first scale J_min_p
  iin = 7;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_axisym_synthesis_mex:InvalidInput:J_min_p",
          "First scale J_min_p must be integer.");
  }
  J_min_p = (int)mxGetScalar(prhs[iin]);
  parameters.J_min_p = J_min_p;

  if (mxGetScalar(prhs[iin]) > (double)J_min_p || J_min_p < 0)
    mexErrMsgIdAndTxt("flaglet_axisym_synthesis_mex:InvalidInput:J_min_p",
          "First scale J_min_p must be positive integer.");

  // Compute ultimate scale J_max
  int J_p = ceil(log(P) / log(B_p));

  if( J_min_p > J_p+1 ) {
    mexErrMsgIdAndTxt("flaglet_axisym_synthesis_mex:InvalidInput:J_min_p",
          "First scale J_min_p must be larger than that!");
  }

  // Parse harmonic band-limit R
  iin = 8;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("slag_synthesis_mex:InvalidInput:Rlimit",
          "Radial limit R must be positive real.");
  }
  parameters.tau = mxGetScalar(prhs[iin]);
  if ( parameters.tau <= 0 )
    mexErrMsgIdAndTxt("slag_synthesis_mex:InvalidInput:RLimitNonInt",
          "Radial limit R must be positive real.");

  // Parse azimuthal band-limit  N
  iin = 9;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_analysis_mex:InvalidInput:N",
          "First scale N must be integer.");
  }
  N = (int)mxGetScalar(prhs[iin]);
  parameters.N = N;

  
  f = (complex double*)calloc( L * (2*L-1) * P, sizeof(complex double));
  flaglet_synthesis(f, f_wav, f_scal, &parameters);
  

  int ntheta = L;
  int nphi = 2 * L - 1;

  // Output function f
  iout = 0;
  plhs[iout] = mxCreateDoubleMatrix(P, L*(2*L-1), mxCOMPLEX);
  f_real = mxGetPr(plhs[iout]);
  f_imag = mxGetPi(plhs[iout]);
  for(n=0; n<P; n++) {    
    for(i=0; i<ntheta*nphi; i++) {
      f_real[i*P+ n] = creal( f[n*ntheta*nphi + i] );
      f_imag[i*P+ n] = cimag( f[n*ntheta*nphi + i] );
    }
  }
  
  free(f);
  free(f_wav);
  free(f_scal);

}
