// FLAGLET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include <flaglet.h>
#include "mex.h"
#include <math.h>
#include <s2let.h>
#include <complex.h>

#define MIN(a,b) ((a) < (b) ? (a) : (b))

/**
 * MATLAB interface: flaglet_axisym_analysis.
 * This function for internal use only.
 * Compute axisymmetric wavelet transform (analysis)
 * with output in pixel space.
 *
 * Usage: 
 *   [f_wav, f_scal] = ...
 *        flaglet_axisym_analysis_mex(f, B_l, B_n, L, N, J_min_l, J_min_n, reality);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  int t, p, n, i, f_m, f_n, reality, downsample;
  int B_l, B_n, L, N, J_min_l, J_min_n;
  double *f_wav_real, *f_scal_real, *f_real, *f_wav_imag, *f_scal_imag, *f_imag;
  complex double *f_wav, *f_scal, *f;
  double *f_wav_r, *f_scal_r, *f_r;
  int iin = 0, iout = 0;

  // Check number of arguments
  if(nrhs!=10) {
    mexErrMsgIdAndTxt("flaglet_axisym_analysis_mex:InvalidInput:nrhs",
          "Require ten inputs.");
  }
  if(nlhs!=2) {
    mexErrMsgIdAndTxt("flaglet_axisym_analysis_mex:InvalidOutput:nlhs",
          "Require two outputs.");
  }

  // Parse reality flag
  iin = 8;
  if( !mxIsLogicalScalar(prhs[iin]) )
    mexErrMsgIdAndTxt("flaglet_axisym_analysis_mex:InvalidInput:reality",
          "Reality flag must be logical.");
  reality = mxIsLogicalScalarTrue(prhs[iin]);

  // Parse multiresolution flag
  iin = 9;
  if( !mxIsLogicalScalar(prhs[iin]) )
    mexErrMsgIdAndTxt("flaglet_axisym_analysis_mex:InvalidInput:downsample",
          "Multiresolution flag must be logical.");
  downsample = mxIsLogicalScalarTrue(prhs[iin]);

  // Parse input dataset f
  iin = 0;
  if( !mxIsDouble(prhs[iin]) ) {
    mexErrMsgIdAndTxt("flag_analysis_mex:InvalidInput:f",
          "Function values must be doubles.");
  }
  f_m = mxGetM(prhs[iin]);
  f_n = mxGetN(prhs[iin]); 
  f_real = mxGetPr(prhs[iin]);
  int f_is_complex = mxIsComplex(prhs[iin]);
  if(reality){
    f_r = (double*)calloc(f_m * f_n, sizeof(double));
    for(t=0; t<f_m; t++)
      for(p=0; p<f_n; p++)
        f_r[t*f_n + p] = f_real[p*f_m + t];
  }
  else {
    f_imag = f_is_complex ? mxGetPi(prhs[iin]) : NULL;
    f = (complex double*)calloc(f_m * f_n, sizeof(complex double));
    for(t=0; t<f_m; t++)
      for(p=0; p<f_n; p++)
        f[t*f_n + p] = f_real[p*f_m + t] 
          + I * (f_is_complex ? f_imag[p*f_m + t] : 0.0);
  }

  // Parse angular wavelet parameter B_l
  iin = 1;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_axisym_analysis_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B_l must be integer.");
  }
  B_l = (int)mxGetScalar(prhs[iin]);

  if (mxGetScalar(prhs[iin]) > (double)B_l || B_l <= 1)
    mexErrMsgIdAndTxt("flaglet_axisym_analysis_mex:InvalidInput:bandLimitNonInt",
          "Wavelet parameter B_l must be positive integer greater than 2");

  // Parse radial wavelet parameter B_n
  iin = 2;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_axisym_analysis_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B_n must be integer.");
  }
  B_n = (int)mxGetScalar(prhs[iin]);

  if (mxGetScalar(prhs[iin]) > (double)B_n || B_n <= 1)
    mexErrMsgIdAndTxt("flaglet_axisym_analysis_mex:InvalidInput:bandLimitNonInt",
          "Wavelet parameter B_n must be positive integer greater than 2");

  // Parse harmonic band-limit L
  iin = 3;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_axisym_analysis_mex:InvalidInput:LbandLimit",
          "Harmonic band-limit L must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);


  if (mxGetScalar(prhs[iin]) > (double)L || L <= 0)
    mexErrMsgIdAndTxt("flaglet_axisym_analysis_mex:InvalidInput:bandLimitNonInt",
          "Harmonic band-limit L must be positive integer.");

  if( B_l >= L ) {
    mexErrMsgIdAndTxt("flaglet_axisym_analysis_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B_l must be smaller than L!");
  }

    // Parse harmonic band-limit N
  iin = 4;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_axisym_analysis_mex:InvalidInput:LbandLimit",
          "Harmonic band-limit N must be integer.");
  }
  N = (int)mxGetScalar(prhs[iin]);


  if (mxGetScalar(prhs[iin]) > (double)N || N <= 0)
    mexErrMsgIdAndTxt("flaglet_axisym_analysis_mex:InvalidInput:bandLimitNonInt",
          "Harmonic band-limit N must be positive integer.");

  if( B_n >= N ) {
    mexErrMsgIdAndTxt("flaglet_axisym_analysis_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B_n must be smaller than N!");
  }
 
  // Parse angular first scale J_min_l
  iin = 5;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_axisym_analysis_mex:InvalidInput:J_min_l",
          "First scale J_min_l must be integer.");
  }
  J_min_l = (int)mxGetScalar(prhs[iin]);

  if (mxGetScalar(prhs[iin]) > (double)J_min_l || J_min_l < 0)
    mexErrMsgIdAndTxt("flaglet_axisym_analysis_mex:InvalidInput:J_min_l",
          "First scale J_min_l must be positive integer.");

  // Compute ultimate scale J_max
  int J_l = ceil(log(L) / log(B_l));


  if( J_min_l > J_l+1 ) {
    mexErrMsgIdAndTxt("flaglet_axisym_analysis_mex:InvalidInput:J_min_l",
          "First scale J_min_l must be larger than that!");
  }

  // Parse angular first scale J_min_n
  iin = 6;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_axisym_analysis_mex:InvalidInput:J_min_n",
          "First scale J_min_n must be integer.");
  }
  J_min_n = (int)mxGetScalar(prhs[iin]);

  if (mxGetScalar(prhs[iin]) > (double)J_min_n || J_min_n < 0)
    mexErrMsgIdAndTxt("flaglet_axisym_analysis_mex:InvalidInput:J_min_n",
          "First scale J_min_n must be positive integer.");

  // Compute ultimate scale J_max
  int J_n = ceil(log(N) / log(B_n));

  if( J_min_n > J_n+1 ) {
    mexErrMsgIdAndTxt("flaglet_axisym_analysis_mex:InvalidInput:J_min_n",
          "First scale J_min_n must be larger than that!");
  }

  // Parse harmonic band-limit R
  iin = 7;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("slag_synthesis_mex:InvalidInput:Rlimit",
          "Radial limit R must be positive real.");
  }
  double R = mxGetScalar(prhs[iin]);
  if ( R <= 0 )
    mexErrMsgIdAndTxt("slag_synthesis_mex:InvalidInput:RLimitNonInt",
          "Radial limit R must be positive real.");

  // Perform wavelet transform in harmonic space and then FLAG reconstruction.
  if(downsample){
    // Multiresolution algorithm
    if(reality){
      flaglet_axisym_allocate_f_wav_multires_real(&f_wav_r, &f_scal_r, B_l, B_n, L, N, J_min_l, J_min_n);
      flaglet_axisym_wav_analysis_multires_real(f_wav_r, f_scal_r, f_r, R, B_l, B_n, L, N, J_min_l, J_min_n);
    }else{
      flaglet_axisym_allocate_f_wav_multires(&f_wav, &f_scal, B_l, B_n, L, N, J_min_l, J_min_n);
      flaglet_axisym_wav_analysis_multires(f_wav, f_scal, f, R, B_l, B_n, L, N, J_min_l, J_min_n); 
    }
  }else{
    // Full-resolution algorithm
    if(reality){
      flaglet_axisym_allocate_f_wav_real(&f_wav_r, &f_scal_r, B_l, B_n, L, N, J_min_l, J_min_n);
      flaglet_axisym_wav_analysis_real(f_wav_r, f_scal_r, f_r, R, B_l, B_n, L, N, J_min_l, J_min_n);
    }else{
      flaglet_axisym_allocate_f_wav(&f_wav, &f_scal, B_l, B_n, L, N, J_min_l, J_min_n);
      flaglet_axisym_wav_analysis(f_wav, f_scal, f, R, B_l, B_n, L, N, J_min_l, J_min_n); 
    }
  }

  // Compute size of wavelet array
  int totalsize = 0;
  if(downsample){
    int jn, jl, bandlimit_n, bandlimit_l;
      for (jn = J_min_n; jn <= J_n; jn++){
        bandlimit_n = MIN(s2let_bandlimit(jn, J_min_n, B_n, N), N);
        for (jl = J_min_l; jl <= J_l; jl++){
          bandlimit_l = MIN(s2let_bandlimit(jl, J_min_l, B_l, L), L);
          totalsize += bandlimit_l * (2 * bandlimit_l - 1) * bandlimit_n;
        }
     }
  }else{
    totalsize = (J_l+1-J_min_l) * L * ( 2 * L - 1 ) * (J_n+1-J_min_n) * N  ;
  }

  // Output wavelets
  if(reality){

    iout = 0;
    plhs[iout] = mxCreateDoubleMatrix(1, totalsize, mxREAL);
    f_wav_real = mxGetPr(plhs[iout]);
    for (i=0; i<totalsize; i++)
      f_wav_real[ i ] = creal( f_wav_r[ i ] );

    iout = 1;
    plhs[iout] = mxCreateDoubleMatrix(N, L*(2*L-1), mxREAL);
    f_scal_real = mxGetPr(plhs[iout]);
    for (n = 0; n < N; n++)
      for (i=0; i<L*(2*L-1); i++)
        f_scal_real[ i * N + n ] = creal(f_scal_r[ n*L*(2*L-1) + i ]);

  }else{

    iout = 0;
    plhs[iout] = mxCreateDoubleMatrix(1, totalsize, mxCOMPLEX);
    f_wav_real = mxGetPr(plhs[iout]);
    f_wav_imag = mxGetPi(plhs[iout]);
    for (i=0; i<totalsize; i++){
      f_wav_real[ i ] = creal( f_wav[ i ] );
      f_wav_imag[ i ] = cimag( f_wav[ i ] );
    }

    iout = 1;
    plhs[iout] = mxCreateDoubleMatrix(N, L*(2*L-1), mxCOMPLEX);
    f_scal_real = mxGetPr(plhs[iout]);
    f_scal_imag = mxGetPi(plhs[iout]);
    for (n = 0; n < N; n++){
      for (i=0; i<L*(2*L-1); i++){
        f_scal_real[ i * N + n ] = creal( f_scal[ n*L*(2*L-1) + i ] );
        f_scal_imag[ i * N + n ] = cimag( f_scal[ n*L*(2*L-1) + i ] );
      }
    }

  }


   if(reality){
    free(f_r);
    free(f_wav_r);
    free(f_scal_r);
  }else{
    free(f);
    free(f_wav);
    free(f_scal);
  }

}
