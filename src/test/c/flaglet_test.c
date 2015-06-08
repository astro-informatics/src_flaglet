// FLAGLET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include "flaglet.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <flag.h>
#include <s2let.h>
#include <ssht.h>

void flaglet_random_flmp(complex double *flmp, int L, int P, int seed)
{
	int i;
	srand( time(NULL) );
	for (i=0; i<L*L*P; i++){
		flmp[i] = (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);
	}
}

void flaglet_random_flmp_real(complex double *flmp, int L, int P, int seed)
{
	int en, el, m, msign, i, i_op, offset;
	int flmsize = L*L;

	for (en=0; en<P; en++) {
		offset = en * flmsize;
		for (el=0; el<L; el++) {
    		m = 0;
    		ssht_sampling_elm2ind(&i, el, m);
    		flmp[offset+i] = (2.0*ran2_dp(seed) - 1.0);
    		for (m=1; m<=el; m++) {
      			ssht_sampling_elm2ind(&i, el, m);
      			flmp[offset+i] = (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);
      			ssht_sampling_elm2ind(&i_op, el, -m);
      			msign = m & 1;
     			msign = 1 - msign - msign; // (-1)^m
     	 		flmp[offset+i_op] = msign * conj(flmp[offset+i]);
    		}
 		}
 	}
}


void flaglet_tiling_test(const flaglet_parameters_t *parameters)
{
	//double kl, kn;
	double *kappa_lp, *kappa0_lp;
	flaglet_tiling_axisym_allocate(&kappa_lp, &kappa0_lp, parameters);

	flaglet_tiling_axisym(kappa_lp, kappa0_lp, parameters);

	/*
	int jl, jp, l, n;
	int J_l = s2let_j_max(L, B_l);
	int J_p = s2let_j_max(P, B_p);
	printf("> KAPPA_0 :       ");
	for (n = 0; n < P; n++){
		for (l = 0; l < L; l++){
			printf(" %f ,", kappa0_lp[n * L + l]);
		}
		printf("\n");}
	printf("\n");
	for (jp = J_min_p; jp <= J_p; jp++){
		for (jl = J_min_l; jl <= J_l; jl++){
			for (n = 0; n < P; n++){
				for (l = 0; l < L; l++){
					printf(" %f ,",kappa_lp[jp*(J_l+1)*L*P  + jl*L*P + n*L + l]);
				}
				printf("\n");}
			printf("\n\n");}
		printf("\n\n\n");}
	*/

	double sum = flaglet_tiling_axisym_check_identity(kappa_lp, kappa0_lp, parameters);
	printf("  - Identity residuals : %10.5e\n", sum);

	free(kappa_lp);
	free(kappa0_lp);
}



void flaglet_lm_test(const flaglet_parameters_t *parameters, int seed)
{
	clock_t time_start, time_end;

	complex double *wav_lmp;
	double *scal_lmp;
	int L = parameters->L;
	int P = parameters->P;

	flaglet_allocate_wav_lmp(&wav_lmp, &scal_lmp, parameters);

	time_start = clock();
	flaglet_wav_lmp(wav_lmp, scal_lmp, parameters);

	printf("  - WAVELET ERROR  : %6.5e\n",
		flaglet_tiling_wavelet_check_identity(wav_lmp, scal_lmp, parameters));

	time_end = clock();
	printf("  - Generate wavelets  : %4.4f seconds\n",
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	complex double *f_wav_lmnp, *f_scal_lmp, *flmp, *flmp_rec;
	flmp = (complex double*)calloc(L * L * P, sizeof(complex double));
	flmp_rec = (complex double*)calloc(L * L * P, sizeof(complex double));

	flaglet_random_flmp(flmp, L, P, seed);

	printf("flaglet_allocate_f_wav_lmnp...");
	flaglet_allocate_f_wav_lmnp(&f_wav_lmnp, &f_scal_lmp, parameters);
	printf("done\n");

	printf("flaglet_analysis_lmnp...");
	time_start = clock();
	flaglet_analysis_lmnp(f_wav_lmnp, f_scal_lmp, flmp, wav_lmp, scal_lmp, parameters);
	time_end = clock();
	printf("  - Wavelet analysis   : %4.4f seconds\n",
		(time_end - time_start) / (double)CLOCKS_PER_SEC);



	time_start = clock();
	flaglet_synthesis_lmnp(flmp_rec, f_wav_lmnp, f_scal_lmp, wav_lmp, scal_lmp, parameters);
	time_end = clock();
	printf("  - Wavelet synthesis  : %4.4f seconds\n",
		(time_end - time_start) / (double)CLOCKS_PER_SEC);


	int n, l, m;
	for (n = 0; n < P; n++){
		for (l = 0; l < L; l++){
	    	for (m = -l; m <= l ; m++){
		    	int ind = n*L*L + l*l + l + m;
		    	if( creal(flmp[ind])-creal(flmp_rec[ind])>0.01 )
					printf("(l,m,p) = (%i,%i,%i) - flmp = (%f,%f) - rec = (%f,%f)\n",l,m,n,creal(flmp[ind]),cimag(flmp[ind]),creal(flmp_rec[ind]),cimag(flmp_rec[ind]));
	   		}
	   	}
	}


	printf("  - Maximum abs error  : %6.5e\n",
		maxerr_cplx(flmp, flmp_rec, L * L * P));

	free(flmp);
	free(flmp_rec);
	free(f_wav_lmnp);
	free(f_scal_lmp);
	free(wav_lmp);
	free(scal_lmp);

}



void flaglet_test(double R, const flaglet_parameters_t *parameters, int seed)
{
	clock_t time_start, time_end;
	int L = parameters->L;
	int P = parameters->P;

	complex double *f, *f_rec, *flmp, *flmp_rec;
	flag_core_allocate_flmn(&flmp, L, P);
	flag_core_allocate_flmn(&flmp_rec, L, P);
	flag_core_allocate_f(&f, L, P);
	flag_core_allocate_f(&f_rec, L, P);

	flaglet_random_flmp(flmp, L, P, seed);

	double *nodes = (double*)calloc(P, sizeof(double));
	double *weights = (double*)calloc(P, sizeof(double));
	assert(nodes != NULL);
	assert(weights != NULL);
	flag_spherlaguerre_sampling(nodes, weights, R, P);

	flag_core_synthesis(f, flmp, nodes, P, L, P);

	free(nodes);
	free(weights);

	complex double *f_wav, *f_scal;
	flaglet_allocate_f_wav(&f_wav, &f_scal, parameters);

	time_start = clock();
	flaglet_analysis(f_wav, f_scal, f, parameters);
	time_end = clock();
	printf("  - Wavelet analysis   : %4.4f seconds\n",
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	time_start = clock();
	flaglet_synthesis(f_rec, f_wav, f_scal, parameters);
	time_end = clock();
	printf("  - Wavelet synthesis  : %4.4f seconds\n",
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	/*printf("  - Maximum abs error  : %6.5e\n",
		maxerr_cplx(f, f_rec, L*(2*L-1)*P));
	for (n = 0; n < L*(2*L-1)*P; n++){
		//printf("f[%i] = (%f,%f) - rec = (%f,%f)\n",n,creal(f[n]),cimag(f[n]),creal(f_rec[n]),cimag(f_rec[n]));
	}*/

	flag_core_analysis(flmp_rec, f_rec, R, L, P);


	int n, l, m;
	for (n = 0; n < P; n++){
		for (l = 0; l < L; l++){
	    	for (m = -l; m <= l ; m++){
		    	int ind = n*L*L + l*l + l + m;
		    	if( creal(flmp[ind])-creal(flmp_rec[ind])>0.01 )
					printf("(l,m,p) = (%i,%i,%i) - flmp = (%f,%f) - rec = (%f,%f)\n",l,m,n,creal(flmp[ind]),cimag(flmp[ind]),creal(flmp_rec[ind]),cimag(flmp_rec[ind]));
	   		}
	   	}
	}

	printf("  - Maximum abs error  : %6.5e\n",
		maxerr_cplx(flmp, flmp_rec, L*L*P));

	/*
	int l, m, n;
	for (n = 0; n < P; n++){
		for (l = 0; l < L; l++){
			for (m = -l; m <= l ; m++){
		    	int ind = n*L*L+l*l+l+m;
				printf("(l,m,n) = (%i,%i,%i) - flmp = (%f,%f) - rec = (%f,%f)\n",l,m,n,creal(flmp[ind]),cimag(flmp[ind]),creal(flmp_rec[ind]),cimag(flmp_rec[ind]));
	    	}
	    }
	}
	*/

	free(f);
	free(f_rec);
	free(f_wav);
	free(f_scal);
}




int main(int argc, char *argv[])
{
	const double R = 10.0;
	const int L = 6;
	const int N = 4;
	const int spin = 0;
	const int P = 6;
	const int B_l = 2;
	const int B_p = 2;
	const int J_min_l = 0;
	const int J_min_p = 0;
	const int upsample = 1;
	const int seed = (int)(10000.0*(double)clock()/(double)CLOCKS_PER_SEC);

	flaglet_parameters_t parameters = {};
	parameters.upsample = upsample;
	parameters.B_l = B_l;
	parameters.B_p = B_p;
	parameters.L = L;
	parameters.P = P;
	parameters.R = R;
	parameters.N = N;
	parameters.spin = spin;
	parameters.J_min_l = J_min_l;
	parameters.J_min_p = J_min_p;

	printf("==========================================================\n");
	printf("PARAMETERS (seed = %i) \n", seed);
	printf(" L = %i  P = %i  N = %i  Bl = %i  Bn = %i  Jminl = %i Jminn = %i \n", L, P, N, B_l, B_p, J_min_l, J_min_p);
	printf("----------------------------------------------------------\n");
	printf("> Testing axisymmetric harmonic tiling...\n");
	flaglet_tiling_test(&parameters);
	printf("==========================================================\n");
	printf("> Testing wavelets in harmonics space...\n");
	flaglet_lm_test(&parameters, seed);
	printf("==========================================================\n");
	printf("> Testing wavelets in pixel space...\n");
	flaglet_test(R, &parameters, seed);
	printf("==========================================================\n");

	return 0;
}
