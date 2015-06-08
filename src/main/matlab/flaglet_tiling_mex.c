// FLAGLET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include <flaglet.h>
#include <s2let.h>
#include "mex.h"

/**
 * MATLAB interface: flaglet_tiling_mex.
 * This function for internal use only.
 * Compute tiling in l-n harmonic space for axisymmetric wavelets.
 *
 * Usage: 
 *   [kappa kappa0] = flaglet_tiling_mex(B_l, B_p, L, P, N, spin, J_min_l, J_min_p);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  int B_l, B_p, L, P, N, spin, J_min_l, J_min_p;
  int iin = 0, iout = 0;

  // Check number of arguments
  if(nrhs!=8) {
    mexErrMsgIdAndTxt("flaglet_tiling_mex:InvalidInput:nrhs",
		      "Require eight inputs.");
  }
  if(nlhs!=2) {
    mexErrMsgIdAndTxt("flaglet_tiling_mex:InvalidOutput:nlhs",
		      "Require two outputs.");
  }

  // Parse angular wavelet parameter B_l
  iin = 0;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_tiling_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B_l must be integer.");
  }
  B_l = (int)mxGetScalar(prhs[iin]);

  if (mxGetScalar(prhs[iin]) > (double)B_l || B_l <= 1)
    mexErrMsgIdAndTxt("flaglet_tiling_mex:InvalidInput:bandLimitNonInt",
          "Wavelet parameter B_l must be positive integer greater than 2");

  // Parse radial wavelet parameter B_p
  iin = 1;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_tiling_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B_p must be integer.");
  }
  B_p = (int)mxGetScalar(prhs[iin]);

  if (mxGetScalar(prhs[iin]) > (double)B_p || B_p <= 1)
    mexErrMsgIdAndTxt("flaglet_tiling_mex:InvalidInput:bandLimitNonInt",
          "Wavelet parameter B_p must be positive integer greater than 2");

  // Parse harmonic band-limit L
  iin = 2;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_tiling_mex:InvalidInput:LbandLimit",
		      "Harmonic band-limit L must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);


  if (mxGetScalar(prhs[iin]) > (double)L || L <= 0)
    mexErrMsgIdAndTxt("flaglet_tiling_mex:InvalidInput:bandLimitNonInt",
		      "Harmonic band-limit L must be positive integer.");

  if( B_l >= L ) {
    mexErrMsgIdAndTxt("flaglet_tiling_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B_l must be smaller than L!");
  }

    // Parse harmonic band-limit P
  iin = 3;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_tiling_mex:InvalidInput:LbandLimit",
          "Harmonic band-limit P must be integer.");
  }
  P = (int)mxGetScalar(prhs[iin]);

  // Parse harmonic band-limit N
  iin = 4;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_tiling_mex:InvalidInput:PbandLimit",
          "Harmonic band-limit N must be integer.");
  }
  N = (int)mxGetScalar(prhs[iin]);

  // Parse harmonic band-limit spin
  iin = 5;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_tiling_mex:InvalidInput:LbandLimit",
          "Harmonic band-limit spin must be integer.");
  }
  spin = (int)mxGetScalar(prhs[iin]);

  // Parse angular first scale J_min_l
  iin = 6;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_tiling_mex:InvalidInput:J_min_l",
          "First scale J_min_l must be integer.");
  }
  J_min_l = (int)mxGetScalar(prhs[iin]);

  if (mxGetScalar(prhs[iin]) > (double)J_min_l || J_min_l < 0)
    mexErrMsgIdAndTxt("flaglet_tiling_mex:InvalidInput:J_min_l",
          "First scale J_min_l must be positive integer.");

  // Compute ultimate scale J_max
  int J_l = flaglet_j_max(L, B_l);

  if( J_min_l > J_l+1 ) {
    mexErrMsgIdAndTxt("flaglet_tiling_mex:InvalidInput:J_min_l",
          "First scale J_min_l must be larger than that!");
  }

  // Parse angular first scale J_min_p
  iin = 7;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flaglet_tiling_mex:InvalidInput:J_min_p",
          "First scale J_min_p must be integer.");
  }
  J_min_p = (int)mxGetScalar(prhs[iin]);

  if (mxGetScalar(prhs[iin]) > (double)J_min_p || J_min_p < 0)
    mexErrMsgIdAndTxt("flaglet_tiling_mex:InvalidInput:J_min_p",
          "First scale J_min_p must be positive integer.");

  // Compute ultimate scale J_max
  int J_p = flaglet_j_max(P, B_p);


  if( J_min_p > J_p+1 ) {
    mexErrMsgIdAndTxt("flaglet_tiling_mex:InvalidInput:J_min_p",
          "First scale J_min_p must be larger than that!");
  }

  flaglet_parameters_t parameters = {};
  parameters.B_l = B_l;
  parameters.B_p = B_p;
  parameters.L = L;
  parameters.P = P;
  parameters.N = N;
  parameters.spin = spin;
  parameters.J_min_l = J_min_l;
  parameters.J_min_p = J_min_p;

  // Allocate arrays
  double *kappa, *kappa0;
  flaglet_tiling_axisym_allocate(&kappa, &kappa0, &parameters);
  flaglet_tiling_axisym(kappa, kappa0, &parameters);

  s2let_parameters_t s2let_parameters = {};
  fill_s2let_angular_parameters(&s2let_parameters, &parameters);

  complex double *s_elm;
  s2let_tiling_direction_allocate(&s_elm, &s2let_parameters);
  s2let_tiling_direction(s_elm, &s2let_parameters);

  // Output kappa and kappa0
  double *kappa_out_real, *kappa_out_imag, *kappa0_out;
  int p, l, jl, jp, m;

  iout = 0;
  plhs[iout] = mxCreateDoubleMatrix(1, (J_l+1)*(J_p+1)*L*L*P, mxCOMPLEX);
  kappa_out_real = mxGetPr(plhs[iout]);
  kappa_out_imag = mxGetPi(plhs[iout]);
  for (jp = J_min_p; jp <= J_p; jp++){
    for (jl = J_min_l; jl <= J_l; jl++){
      for(p=0; p<P; p++){
        for (l = 0; l < L; l++){
          //if (kappa[ jp*(J_l+1)*L*P  + jl*L*P + p*L + l ] > 0)
          //  printf("ind, jp,jl,p,l = %i, %i,%i,%i,%i   kappa=(%f)\n", jp*(J_l+1)*L*L*P  + jl*L*L*P + p*L*L + l*l + l + m, jp,jl,p,l, kappa[ jp*(J_l+1)*L*P  + jl*L*P + p*L + l ]);

          for (m = -l; m <= l; m++){
            //printf("ind, jp,jl,p,l,m = %i, %i,%i,%i,%i,%i   kappa=(%f)   selm=(%f,%f)\n", jp*(J_l+1)*L*L*P  + jl*L*L*P + p*L*L + l*l + l + m, jp,jl,p,l,m, kappa[ jp*(J_l+1)*L*P  + jl*L*P + p*L + l ], creal(s_elm[l*l + l + m]), cimag(s_elm[l*l + l + m]));
            kappa_out_real[ jp*(J_l+1)*L*L*P  + jl*L*L*P + p*L*L + l*l + l + m ] = 
              kappa[ jp*(J_l+1)*L*P  + jl*L*P + p*L + l ] * creal(s_elm[l*l + l + m]);
            kappa_out_imag[ jp*(J_l+1)*L*L*P  + jl*L*L*P + p*L*L + l*l + l + m ] = 
              kappa[ jp*(J_l+1)*L*P  + jl*L*P + p*L + l ] * cimag(s_elm[l*l + l + m]);
            }
          }
        }
      }
    }

  iout = 1;
  plhs[iout] = mxCreateDoubleMatrix(1, L*L*P, mxREAL);
  kappa0_out = mxGetPr(plhs[iout]);
  for(p=0; p<P; p++){
    for(l=0; l<L; l++){
      for(m=-l; m<=l; m++){
        //printf("ind,ind2,p,l,m = %i,%i,%i,%i,%i  kappa0=%f\n",p*L*L + l*l + l + m, p * L + l ,p,l,m,kappa0[ p * L + l ]);
        kappa0_out[ p*L*L + l*l + l + m ] = creal(kappa0[ p * L + l ]); 
      }
    }
  }

  free(kappa);
  free(kappa0);

}
