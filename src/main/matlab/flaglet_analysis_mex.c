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
 * MATLAB interface: flaglet_analysis.
 * This function for internal use only.
 * Compute axisymmetric wavelet transform (analysis)
 * with output in pixel space.
 *
 * Usage: 
 *   [f_wav, f_scal] = ...
 *        flaglet_analysis_mex(f, B_l, B_p, L, P, J_min_l, J_min_p, N, downsample);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  int t, p, n, i, f_m, f_n;
  s2let_parameters_t parameters_ang = {};
  s2let_parameters_t parameters_rad = {};
  flaglet_parameters_t parameters = {};
  double *f_wav_real, *f_scal_real, *f_real, *f_wav_imag, *f_scal_imag, *f_imag;
  int B_l, B_p, L, P, Nj, N, J_min_l, J_min_p;
  complex double *f_wav, *f_scal, *f;
  int iin = 0, iout = 0;

  // Check number of arguments
  if(nrhs!=10) {
    mexErrMsgIdAndTxt("flaglet_analysis_mex:InvalidInput:nrhs",
          "Require ten inputs.");
  }
  if(nlhs!=2) {
    mexErrMsgIdAndTxt("flaglet_analysis_mex:InvalidOutput:nlhs",
          "Require two outputs.");
  }

  // Parse multiresolution flag
  iin = 9;
  if( !mxIsLogicalScalar(prhs[iin]) )
    mexErrMsgIdAndTxt("flaglet_analysis_mex:InvalidInput:upsample",
          "Multiresolution flag must be logical.");
  parameters.upsample = !mxIsLogicalScalarTrue(prhs[iin]);

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
  f_imag = mxGetPi(prhs[iin]);
  f = (complex double*)calloc(f_m * f_n, sizeof(complex double));
  for(t=0; t<f_m; t++)
    for(p=0; p<f_n; p++)
      f[t*f_n + p] = f_real[p*f_m + t] 
        + I * (f_is_complex ? f_imag[p*f_m + t] : 0.0);

  // Parse angular wavelet parameter B_l
  iin = 1;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_analysis_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B_l must be integer.");
  }
  B_l = (int)mxGetScalar(prhs[iin]);
  parameters.B_l = B_l;

  if (mxGetScalar(prhs[iin]) > (double)B_l || B_l <= 1)
    mexErrMsgIdAndTxt("flaglet_analysis_mex:InvalidInput:bandLimitNonInt",
          "Wavelet parameter B_l must be positive integer greater than 2");

  // Parse radial wavelet parameter B_p
  iin = 2;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_analysis_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B_p must be integer.");
  }
  B_p = (int)mxGetScalar(prhs[iin]);
  parameters.B_p = B_p;

  if (mxGetScalar(prhs[iin]) > (double)B_p || B_p <= 1)
    mexErrMsgIdAndTxt("flaglet_analysis_mex:InvalidInput:bandLimitNonInt",
          "Wavelet parameter B_p must be positive integer greater than 2");

  // Parse harmonic band-limit L
  iin = 3;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_analysis_mex:InvalidInput:LbandLimit",
          "Harmonic band-limit L must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);
  parameters.L = L;

  if (mxGetScalar(prhs[iin]) > (double)L || L <= 0)
    mexErrMsgIdAndTxt("flaglet_analysis_mex:InvalidInput:bandLimitNonInt",
          "Harmonic band-limit L must be positive integer.");

  if( B_l >= L ) {
    mexErrMsgIdAndTxt("flaglet_analysis_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B_l must be smaller than L!");
  }

    // Parse harmonic band-limit P
  iin = 4;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_analysis_mex:InvalidInput:LbandLimit",
          "Harmonic band-limit P must be integer.");
  }
  P = (int)mxGetScalar(prhs[iin]);
  parameters.P = P;

  if (mxGetScalar(prhs[iin]) > (double)P || P <= 0)
    mexErrMsgIdAndTxt("flaglet_analysis_mex:InvalidInput:bandLimitNonInt",
          "Harmonic band-limit P must be positive integer.");

  if( B_p >= P ) {
    mexErrMsgIdAndTxt("flaglet_analysis_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B_p must be smaller than P!");
  }
 
  // Parse angular first scale J_min_l
  iin = 5;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_analysis_mex:InvalidInput:J_min_l",
          "First scale J_min_l must be integer.");
  }
  J_min_l = (int)mxGetScalar(prhs[iin]);
  parameters.J_min_l = J_min_l;

  if (mxGetScalar(prhs[iin]) > (double)J_min_l || J_min_l < 0)
    mexErrMsgIdAndTxt("flaglet_analysis_mex:InvalidInput:J_min_l",
          "First scale J_min_l must be positive integer.");

  // Compute ultimate scale J_max
  int J_l = ceil(log(L) / log(B_l));


  if( J_min_l > J_l+1 ) {
    mexErrMsgIdAndTxt("flaglet_analysis_mex:InvalidInput:J_min_l",
          "First scale J_min_l must be larger than that!");
  } 

  // Parse angular first scale J_min_p
  iin = 6;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_analysis_mex:InvalidInput:J_min_p",
          "First scale J_min_p must be integer.");
  }
  J_min_p = (int)mxGetScalar(prhs[iin]);
  parameters.J_min_p = J_min_p;

  if (mxGetScalar(prhs[iin]) > (double)J_min_p || J_min_p < 0)
    mexErrMsgIdAndTxt("flaglet_analysis_mex:InvalidInput:J_min_p",
          "First scale J_min_p must be positive integer.");

  // Compute ultimate scale J_max
  int J_p = ceil(log(P) / log(B_p));

  if( J_min_p > J_p+1 ) {
    mexErrMsgIdAndTxt("flaglet_analysis_mex:InvalidInput:J_min_p",
          "First scale J_min_p must be larger than that!");
  }

  // Parse harmonic band-limit R
  iin = 7;
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
  iin = 8;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_analysis_mex:InvalidInput:N",
          "First scale N must be integer.");
  }
  N = (int)mxGetScalar(prhs[iin]);
  parameters.N = N;


  parameters_ang.L = L;
  parameters_ang.J_min = J_min_l;
  parameters_ang.B = B_l;
  parameters_rad.L = P;
  parameters_rad.J_min = J_min_p;
  parameters_rad.B = B_p;

  so3_parameters_t so3_parameters = {};
  fill_so3_angular_parameters(&so3_parameters, &parameters);

  // Compute size of wavelet array
  int totalsize = 0;
  if(!parameters.upsample){
    int jn, jl, bandlimit_p, bandlimit_l;
      for (jn = J_min_p; jn <= J_p; jn++){
        bandlimit_p = MIN(s2let_bandlimit(jn, &parameters_rad), P);
        for (jl = J_min_l; jl <= J_l; jl++){
          bandlimit_l = MIN(s2let_bandlimit(jl, &parameters_ang), L);
          so3_parameters.L = bandlimit_l;
          so3_parameters.N = N;
          totalsize += bandlimit_p * so3_sampling_f_size(&so3_parameters);
        }
     }
  }else{
    totalsize = (J_l+1-J_min_l) * (J_p+1-J_min_p) * L * (2*L-1) * P * N  ;
  }

  // Perform wavelet transform in harmonic space and then FLAG reconstruction.
  //f_wav = (complex double*)calloc( totalsize, sizeof(complex double));
  //f_scal = (complex double*)calloc( L * (2*L-1) * P, sizeof(complex double));
  flaglet_allocate_f_wav(&f_wav, &f_scal, &parameters);
  flaglet_analysis(f_wav, f_scal, f, &parameters);

  iout = 0;
  plhs[iout] = mxCreateDoubleMatrix(1, totalsize, mxCOMPLEX);
  f_wav_real = mxGetPr(plhs[iout]);
  f_wav_imag = mxGetPi(plhs[iout]);
  for (i=0; i<totalsize; i++){
    f_wav_real[ i ] = creal( f_wav[ i ] );
    f_wav_imag[ i ] = cimag( f_wav[ i ] );
  }

  iout = 1;
  plhs[iout] = mxCreateDoubleMatrix(P, L*(2*L-1), mxCOMPLEX);
  f_scal_real = mxGetPr(plhs[iout]);
  f_scal_imag = mxGetPi(plhs[iout]);
  for (n = 0; n < P; n++){
    for (i=0; i<L*(2*L-1); i++){
      f_scal_real[ i * P + n ] = creal( f_scal[ n*L*(2*L-1) + i ] );
      f_scal_imag[ i * P + n ] = cimag( f_scal[ n*L*(2*L-1) + i ] );
    }
  }

  free(f);
  free(f_wav);
  free(f_scal);

}
