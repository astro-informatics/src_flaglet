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
#include <flag/flag.h>
#include <s2let/s2let.h>
#include <ssht/ssht.h>

double dot_test(complex double *a, complex double *b, int size)
{
	double value = 0;
	int i;
	for(i = 0; i<size; i++){
		value += conj(a[i])*b[i];
	}
	return cabs(value);
}

double dot_test_real(double *a, double *b, int size)
{
	double value = 0;
	int i;
	for(i = 0; i<size; i++){
		value += a[i]*b[i];
	}
	return value;
}

void flaglet_random_flmp(complex double *flmp, int L, int P, int spin, int seed)
{
	int i, index, en, el, m, offset;
	srand( time(NULL) );
	for (i=0; i<L*L*P; i++){
		flmp[i] = (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);
	}
	for (en=0; en<P; en++) {
		offset = en * L * L;
		for (el=0; el<=abs(spin); el++) {
    		for (m=-el; m<=el; m++) {
    			ssht_sampling_elm2ind(&index, el, m);
    			flmp[offset+index] = 0;
    		}
    	}
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
	int spin = parameters->spin;

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

	flaglet_random_flmp(flmp, L, P, spin, seed);

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



void flaglet_test(double tau, const flaglet_parameters_t *parameters, int seed)
{
	clock_t time_start, time_end;
	int L = parameters->L;
	int P = parameters->P;
	int spin = parameters->spin;

	complex double *f, *f_rec, *flmp, *flmp_rec;
	flag_core_allocate_flmn(&flmp, L, P);
	flag_core_allocate_flmn(&flmp_rec, L, P);
	flag_core_allocate_f(&f, L, P);
	flag_core_allocate_f(&f_rec, L, P);

	flaglet_random_flmp(flmp, L, P, spin, seed);

	double *nodes = (double*)calloc(P, sizeof(double));
	double *weights = (double*)calloc(P, sizeof(double));
	assert(nodes != NULL);
	assert(weights != NULL);
	flag_spherlaguerre_sampling(nodes, weights, tau, P);

	flag_core_synthesis(f, flmp, nodes, P, L, tau, P, spin);

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

	flag_core_analysis(flmp_rec, f_rec, L, tau, P, spin);

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

	free(f);
	free(f_rec);
	free(f_wav);
	free(f_scal);
}


void flaglet_analysis_adjoint_test(double tau, const flaglet_parameters_t *parameters, int seed)
{
	clock_t time_start, time_end;
	int L = parameters->L;
	int P = parameters->P;
	int spin = parameters->spin;

	complex double *f_input, *f_output, *f, *flmp_temp, *flmp_input, *flmp_output;
	flag_core_allocate_flmn(&flmp_temp, L, P);
	flag_core_allocate_flmn(&flmp_input, L, P);
	flag_core_allocate_f(&f, L, P);
	flag_core_allocate_f(&f_input, L, P);
	flag_core_allocate_f(&f_output, L, P);

	flaglet_random_flmp(flmp_temp, L, P, spin, seed);
	flaglet_random_flmp(flmp_input, L, P, spin, seed);

	double *nodes = (double*)calloc(P, sizeof(double));
	double *weights = (double*)calloc(P, sizeof(double));
	assert(nodes != NULL);
	assert(weights != NULL);
	flag_spherlaguerre_sampling(nodes, weights, tau, P);
	flag_core_synthesis(f_input, flmp_input, nodes, P, L, tau, P, spin);
	flag_core_synthesis(f, flmp_temp, nodes, P, L, tau, P, spin);

	free(nodes);
	free(weights);
	free(flmp_input);
	free(flmp_temp);

	complex double *f_wav_input, *f_scal_input;
	flaglet_allocate_f_wav(&f_wav_input, &f_scal_input, parameters);
	flaglet_analysis(f_wav_input, f_scal_input, f, parameters);
	flaglet_analysis_adjoint(f_output, f_wav_input, f_scal_input, parameters);

	complex double *f_wav_output, *f_scal_output;
	flaglet_allocate_f_wav(&f_wav_output, &f_scal_output, parameters);
	flaglet_analysis(f_wav_output, f_scal_output, f_input, parameters);

	printf("  - Dot error on ADJOINT ANALYSIS FLAGLET operators  : %6.5e\n",
		dot_test(f_wav_input, f_wav_output, flaglet_wav_size(parameters))
	  + dot_test(f_scal_input, f_scal_output, flaglet_scal_size(parameters))
	  - dot_test(f_input, f_output, flaglet_f_size(parameters)));

	free(f_input);
	free(f_output);
	free(f);
	free(f_wav_input);
	free(f_scal_input);
	free(f_wav_output);
	free(f_scal_output);
}

void flaglet_synthesis_adjoint_test(double tau, const flaglet_parameters_t *parameters, int seed)
{
	clock_t time_start, time_end;
	int L = parameters->L;
	int P = parameters->P;
	int spin = parameters->spin;

	complex double *f_input, *f_output, *f, *flmp_temp, *flmp_input, *flmp_output;
	flag_core_allocate_flmn(&flmp_temp, L, P);
	flag_core_allocate_flmn(&flmp_input, L, P);
	flag_core_allocate_f(&f, L, P);
	flag_core_allocate_f(&f_input, L, P);
	flag_core_allocate_f(&f_output, L, P);

	flaglet_random_flmp(flmp_temp, L, P, spin, seed);
	flaglet_random_flmp(flmp_input, L, P, spin, seed);

	double *nodes = (double*)calloc(P, sizeof(double));
	double *weights = (double*)calloc(P, sizeof(double));
	assert(nodes != NULL);
	assert(weights != NULL);
	flag_spherlaguerre_sampling(nodes, weights, tau, P);
	flag_core_synthesis(f_input, flmp_input, nodes, P, L, tau, P, spin);
	flag_core_synthesis(f, flmp_temp, nodes, P, L, tau, P, spin);

	free(nodes);
	free(weights);
	free(flmp_input);
	free(flmp_temp);

	complex double *f_wav_input, *f_scal_input;
	flaglet_allocate_f_wav(&f_wav_input, &f_scal_input, parameters);
	flaglet_analysis(f_wav_input, f_scal_input, f, parameters);
	flaglet_synthesis(f_output, f_wav_input, f_scal_input, parameters);

	complex double *f_wav_output, *f_scal_output;
	flaglet_allocate_f_wav(&f_wav_output, &f_scal_output, parameters);
	flaglet_synthesis_adjoint(f_wav_output, f_scal_output, f_input, parameters);

	printf("  - Dot error on ADJOINT SYNTHESIS FLAGLET operators  : %6.5e\n",
		dot_test(f_wav_input, f_wav_output, flaglet_wav_size(parameters))
	  + dot_test(f_scal_input, f_scal_output, flaglet_scal_size(parameters))
	  - dot_test(f_input, f_output, flaglet_f_size(parameters)));

	free(f_input);
	free(f_output);
	free(f);
	free(f_wav_input);
	free(f_scal_input);
	free(f_wav_output);
	free(f_scal_output);
}




int main(int argc, char *argv[])
{
	const double tau = 1.0;
	const int L = 32;
	const int N = 4;
	const int spin = 2;
	const int P = 32;
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
	parameters.tau = tau;
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
	printf("----------------------------------------------------------\n");
	printf("> Testing wavelets in harmonics space...\n");
	flaglet_lm_test(&parameters, seed);
	printf("----------------------------------------------------------\n");
	printf("> Testing wavelets in pixel space...\n");
	flaglet_test(tau, &parameters, seed);
	printf("==========================================================\n");
	printf("> Testing wavelets analysis adjoints in pixel space...\n");
	flaglet_analysis_adjoint_test(tau, &parameters, seed);
	printf("==========================================================\n");
	printf("> Testing wavelets synthesis adjoints in pixel space...\n");
	flaglet_synthesis_adjoint_test(tau, &parameters, seed);
	printf("==========================================================\n");

	return 0;
}
