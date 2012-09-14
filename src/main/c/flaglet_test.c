// FLAGLET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "flaglet.h"
#include <assert.h>

#define MAX(a,b) ((a) > (b) ? (a) : (b))

double maxerr_cplx(complex double *a, complex double *b, int size)
{
	double value = 0;
	int i;
	for(i = 0; i<size; i++){
		value = MAX( cabs( a[i]-b[i] ), value );
	}
	return value;
}

double ran2_dp(int idum) {
  int IM1=2147483563,IM2=2147483399,IMM1=IM1-1, 
    IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, 
    PTAB=32,PDIV=1+IMM1/PTAB;

  double AM=1./IM1,EPS=1.2e-7,RPMX=1.-EPS;
  int j,k;
  static int iv[32],iy,idum2 = 123456789; 
  // P.B. in C static variables are initialised to 0 by default.

  if (idum <= 0) {
    idum= (-idum>1 ? -idum : 1); // max(-idum,1);
    idum2=idum;
    for(j=PTAB+8;j>=1;j--) {
      k=idum/IQ1;
      idum=IA1*(idum-k*IQ1)-k*IR1;
      if (idum < 0) idum=idum+IM1;
      if (j < PTAB) iv[j-1]=idum;
    }
    iy=iv[0];
  }
  k=idum/IQ1;
  idum=IA1*(idum-k*IQ1)-k*IR1;
  if (idum < 0) idum=idum+IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2=idum2+IM2;
  j=1+iy/PDIV;
  iy=iv[j-1]-idum2;
  iv[j-1]=idum;
  if(iy < 1)iy=iy+IMM1;
  return (AM*iy < RPMX ? AM*iy : RPMX); // min(AM*iy,RPMX);
}

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
	int flmsize = ssht_flm_size(L);

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

void flaglet_tilling_test(int B_l, int B_p, int L, int P, int J_min_l, int J_min_p)
{		
	//double kl, kn;
	double *kappa_lp, *kappa0_lp;
	flaglet_axisym_allocate_tilling(&kappa_lp, &kappa0_lp, B_l, B_p, L, P);
	
	flaglet_axisym_tilling(kappa_lp, kappa0_lp, B_l, B_p, L, P, J_min_l, J_min_p);
	
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

	double sum = flaglet_axisym_check_identity(kappa_lp, kappa0_lp, B_l, B_p, L, P, J_min_l, J_min_p);
	printf("  - Identity residuals : %e\n", sum);

	free(kappa_lp);
	free(kappa0_lp);
}

void flaglet_axisym_wav_lm_test(int B_l, int B_p, int L, int P, int J_min_l, int J_min_p, int seed)
{		
	clock_t time_start, time_end;

	double *wav_lmp, *scal_lmp;

	flaglet_axisym_allocate_wav_lmp(&wav_lmp, &scal_lmp, B_l, B_p, L, P);

	time_start = clock();
	flaglet_axisym_wav_lmp(wav_lmp, scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);
	time_end = clock();
	printf("  - Generate wavelets  : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	complex double *f_wav_lmp, *f_scal_lmp, *flmp, *flmp_rec;
	flmp = (complex double*)calloc(L * L * P, sizeof(complex double));
	flmp_rec = (complex double*)calloc(L * L * P, sizeof(complex double));

	flaglet_random_flmp(flmp, L, P, seed);
	
	flaglet_axisym_allocate_f_wav_lmp(&f_wav_lmp, &f_scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);

	time_start = clock();
	flaglet_axisym_wav_analysis_lmp(f_wav_lmp, f_scal_lmp, flmp, wav_lmp, scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);
	time_end = clock();
	printf("  - Wavelet analysis   : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	/*
	int n, l, m;
	for (n = 0; n < P; n++){
		for (l = 0; l < L; l++){
	    	for (m = -l; m <= l ; m++){
		    	int ind = n*L*L + l*l + l + m;
				printf("(l,m,n) = (%i,%i,%i) - f_scal_lmp = (%f,%f)\n",l,m,n,creal(f_scal_lmp[ind]),cimag(f_scal_lmp[ind]));
	   		}
	   	}
	}
	*/

	time_start = clock();
	flaglet_axisym_wav_synthesis_lmp(flmp_rec, f_wav_lmp, f_scal_lmp, wav_lmp, scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);
	time_end = clock();
	printf("  - Wavelet synthesis  : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	/*
	for (n = 0; n < P; n++){
		for (l = 0; l < L; l++){
	    	for (m = -l; m <= l ; m++){
		    	int ind = n*L*L + l*l + l + m;
		    	if( creal(flmp[ind])-creal(flmp_rec[ind])>0.01 )
					printf("(l,m,n) = (%i,%i,%i) - flmp = (%f,%f) - rec = (%f,%f)\n",l,m,n,creal(flmp[ind]),cimag(flmp[ind]),creal(flmp_rec[ind]),cimag(flmp_rec[ind]));
	   		}
	   	}
	}
	*/

	printf("  - Maximum abs error  : %6.5e\n", 
		maxerr_cplx(flmp, flmp_rec, L * L * P));

	free(flmp);
	free(flmp_rec);
	free(f_wav_lmp);
	free(f_scal_lmp);
	free(wav_lmp);
	free(scal_lmp);

}


void flaglet_axisym_wav_lm_multires_test(int B_l, int B_p, int L, int P, int J_min_l, int J_min_p, int seed)
{		
	clock_t time_start, time_end;

	double *wav_lmp, *scal_lmp;

	flaglet_axisym_allocate_wav_lmp(&wav_lmp, &scal_lmp, B_l, B_p, L, P);

	time_start = clock();
	flaglet_axisym_wav_lmp(wav_lmp, scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);
	time_end = clock();
	printf("  - Generate wavelets  : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	complex double *f_wav_lmp, *f_scal_lmp, *flmp, *flmp_rec;
	flmp = (complex double*)calloc(L * L * P, sizeof(complex double));
	flmp_rec = (complex double*)calloc(L * L * P, sizeof(complex double));

	flaglet_random_flmp(flmp, L, P, seed);
	
	flaglet_axisym_allocate_f_wav_multires_lmp(&f_wav_lmp, &f_scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);

	time_start = clock();
	flaglet_axisym_wav_analysis_multires_lmp(f_wav_lmp, f_scal_lmp, flmp, wav_lmp, scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);
	time_end = clock();
	printf("  - Wavelet analysis   : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	time_start = clock();
	flaglet_axisym_wav_synthesis_multires_lmp(flmp_rec, f_wav_lmp, f_scal_lmp, wav_lmp, scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);
	time_end = clock();
	printf("  - Wavelet synthesis  : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	printf("  - Maximum abs error  : %6.5e\n", 
		maxerr_cplx(flmp, flmp_rec, L * L * P));

	/*
	int n, l, m;
	for (n = 0; n < P; n++){
		for (l = 0; l < L; l++){
	    	for (m = -l; m <= l ; m++){
		    	int ind = n*L*L + l*l + l + m;
		    	if( creal(flmp[ind])-creal(flmp_rec[ind])>0.01 )
					printf("(l,m,n) = (%i,%i,%i) - flmp = (%f,%f) - rec = (%f,%f)\n",l,m,n,creal(flmp[ind]),cimag(flmp[ind]),creal(flmp_rec[ind]),cimag(flmp_rec[ind]));
	   		}
	   	}
	}
	*/

	free(flmp);
	free(flmp_rec);
	free(f_wav_lmp);
	free(f_scal_lmp);
	free(wav_lmp);
	free(scal_lmp);

}


void flaglet_axisym_wav_test(double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p, int seed)
{
	clock_t time_start, time_end;

	complex double *f, *f_rec, *flmp, *flmp_rec;
	flag_allocate_flmn(&flmp, L, P);
	flag_allocate_flmn(&flmp_rec, L, P);
	flag_allocate_f(&f, L, P);
	flag_allocate_f(&f_rec, L, P);

	flaglet_random_flmp(flmp, L, P, seed);

	double *nodes = (double*)calloc(P, sizeof(double));
	double *weights = (double*)calloc(P, sizeof(double));
	assert(nodes != NULL);
	assert(weights != NULL);
	flag_spherlaguerre_sampling(nodes, weights, R, P);

	flag_synthesis(f, flmp, nodes, P, L, P);

	free(nodes);
	free(weights);

	complex double *f_wav, *f_scal;
	flaglet_axisym_allocate_f_wav(&f_wav, &f_scal, B_l, B_p, L, P, J_min_l, J_min_p);

	time_start = clock();
	flaglet_axisym_wav_analysis(f_wav, f_scal, f, R, B_l, B_p, L, P, J_min_l, J_min_p);
	time_end = clock();
	printf("  - Wavelet analysis   : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	time_start = clock();
	flaglet_axisym_wav_synthesis(f_rec, f_wav, f_scal, R, B_l, B_p, L, P, J_min_l, J_min_p);
	time_end = clock();
	printf("  - Wavelet synthesis  : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	/*printf("  - Maximum abs error  : %6.5e\n", 
		maxerr_cplx(f, f_rec, L*(2*L-1)*P));
	for (n = 0; n < L*(2*L-1)*P; n++){
		//printf("f[%i] = (%f,%f) - rec = (%f,%f)\n",n,creal(f[n]),cimag(f[n]),creal(f_rec[n]),cimag(f_rec[n]));
	}*/

	flag_analysis(flmp_rec, f_rec, R, L, P);

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

void flaglet_axisym_wav_real_test(double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p, int seed)
{
	clock_t time_start, time_end;

	complex *flmp, *flmp_rec;
	double *f, *f_rec;
	flag_allocate_flmn(&flmp, L, P);
	flag_allocate_flmn(&flmp_rec, L, P);
	flag_allocate_f_real(&f, L, P);
	flag_allocate_f_real(&f_rec, L, P);

	flaglet_random_flmp_real(flmp, L, P, seed);

	double *nodes = (double*)calloc(P, sizeof(double));
	double *weights = (double*)calloc(P, sizeof(double));
	assert(nodes != NULL);
	assert(weights != NULL);
	flag_spherlaguerre_sampling(nodes, weights, R, P);

	flag_synthesis_real(f, flmp, nodes, P, L, P);

	free(nodes);
	free(weights);

	double *f_wav, *f_scal;
	flaglet_axisym_allocate_f_wav_real(&f_wav, &f_scal, B_l, B_p, L, P, J_min_l, J_min_p);

	time_start = clock();
	flaglet_axisym_wav_analysis_real(f_wav, f_scal, f, R, B_l, B_p, L, P, J_min_l, J_min_p);
	time_end = clock();
	printf("  - Wavelet analysis   : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	time_start = clock();
	flaglet_axisym_wav_synthesis_real(f_rec, f_wav, f_scal, R, B_l, B_p, L, P, J_min_l, J_min_p);
	time_end = clock();
	printf("  - Wavelet synthesis  : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	/*
	printf("  - Maximum abs error  : %6.5e\n", 
		maxerr_cplx(f, f_rec, L*(2*L-1)*P));
	for (n = 0; n < L*(2*L-1)*P; n++){
		printf("f[%i] = (%f,%f) - rec = (%f,%f)\n",n,creal(f[n]),cimag(f[n]),creal(f_rec[n]),cimag(f_rec[n]));
	}
	*/

	flag_analysis_real(flmp_rec, f_rec, R, L, P);

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


void flaglet_axisym_wav_multires_test(double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p, int seed)
{
	clock_t time_start, time_end;

	complex double *f, *f_rec, *flmp, *flmp_rec;
	flag_allocate_flmn(&flmp, L, P);
	flag_allocate_flmn(&flmp_rec, L, P);
	flag_allocate_f(&f, L, P);
	flag_allocate_f(&f_rec, L, P);

	flaglet_random_flmp(flmp, L, P, seed);

	double *nodes = (double*)calloc(P, sizeof(double));
	double *weights = (double*)calloc(P, sizeof(double));
	assert(nodes != NULL);
	assert(weights != NULL);
	flag_spherlaguerre_sampling(nodes, weights, R, P);

	flag_synthesis(f, flmp, nodes, P, L, P);

	free(nodes);
	free(weights);

	complex double *f_wav, *f_scal;
	flaglet_axisym_allocate_f_wav_multires(&f_wav, &f_scal, B_l, B_p, L, P, J_min_l, J_min_p);

	time_start = clock();
	flaglet_axisym_wav_analysis_multires(f_wav, f_scal, f, R, B_l, B_p, L, P, J_min_l, J_min_p);
	time_end = clock();
	printf("  - Wavelet analysis   : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	time_start = clock();
	flaglet_axisym_wav_synthesis_multires(f_rec, f_wav, f_scal, R, B_l, B_p, L, P, J_min_l, J_min_p);
	time_end = clock();
	printf("  - Wavelet synthesis  : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	flag_analysis(flmp_rec, f_rec, R, L, P);

	printf("  - Maximum abs error  : %6.5e\n", 
		maxerr_cplx(flmp, flmp_rec, L*L*P));
	
	/*
	int l, m, n;
	for (n = 0; n < P; n++){
		for (l = 0; l < L; l++){
			for (m = -l; m <= l ; m++){
		    	int ind = n*L*L+l*l+l+m;
		    	if( creal(flmp[ind])-creal(flmp_rec[ind])>0.01 )
					printf("(l,m,n) = (%i,%i,%i) - flmp = (%f,%f) - rec = (%f,%f)\n",l,m,n,creal(flmp[ind]),cimag(flmp[ind]),creal(flmp_rec[ind]),cimag(flmp_rec[ind]));
	    		
	    	}
	    	//printf(" %f ", flmp[n*L*L+l*l+l] - flmp_rec[n*L*L+l*l+l]);
	    }
	    printf("\n");
	}
	*/
	
	free(f);
	free(f_rec);
	free(f_wav);
	free(f_scal);
}

void flaglet_axisym_wav_multires_real_test(double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p, int seed)
{
	clock_t time_start, time_end;
	int J_l = s2let_j_max(L, B_l);
	int J_p = s2let_j_max(P, B_p);

	complex *flmp, *flmp_rec;
	double *f, *f_rec;
	flag_allocate_flmn(&flmp, L, P);
	flag_allocate_flmn(&flmp_rec, L, P);
	flag_allocate_f_real(&f, L, P);
	flag_allocate_f_real(&f_rec, L, P);

	flaglet_random_flmp_real(flmp, L, P, seed);

	double *nodes = (double*)calloc(P, sizeof(double));
	double *weights = (double*)calloc(P, sizeof(double));
	assert(nodes != NULL);
	assert(weights != NULL);
	flag_spherlaguerre_sampling(nodes, weights, R, P);

	flag_synthesis_real(f, flmp, nodes, P, L, P);

	free(nodes);
	free(weights);

	double *f_wav, *f_scal;
	flaglet_axisym_allocate_f_wav_multires_real(&f_wav, &f_scal, B_l, B_p, L, P, J_min_l, J_min_p);

	time_start = clock();
	flaglet_axisym_wav_analysis_multires_real(f_wav, f_scal, f, R, B_l, B_p, L, P, J_min_l, J_min_p);
	time_end = clock();
	printf("  - Wavelet analysis   : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	time_start = clock();
	flaglet_axisym_wav_synthesis_multires_real(f_rec, f_wav, f_scal, R, B_l, B_p, L, P, J_min_l, J_min_p);
	time_end = clock();
	printf("  - Wavelet synthesis  : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	/*
	printf("  - Maximum abs error  : %6.5e\n", 
		maxerr_cplx(f, f_rec, L*(2*L-1)*P));
	for (n = 0; n < L*(2*L-1)*P; n++){
		printf("f[%i] = (%f,%f) - rec = (%f,%f)\n",n,creal(f[n]),cimag(f[n]),creal(f_rec[n]),cimag(f_rec[n]));
	}
	*/

	flag_analysis_real(flmp_rec, f_rec, R, L, P);

	printf("  - Maximum abs error  : %6.5e\n", 
		maxerr_cplx(flmp, flmp_rec, L*L*P));
	
	/*
	int l, m, n;
	for (n = 0; n < P; n++){
		for (l = 0; l < L; l++){
			for (m = -l; m <= l ; m++){
		    	int ind = n*L*L+l*l+l+m;
		    	if( creal(flmp[ind])-creal(flmp_rec[ind])>0.001 )
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


void flaglet_axisym_wav_performance_test(double R, int NREPEAT, int NSCALE, int seed, int B_l, int B_p, int J_min_l, int J_min_p)
{
	complex double *f, *f_rec, *flmp, *flmp_rec;
	clock_t time_start, time_end;
	int sc, repeat;
	double tottime_analysis = 0, tottime_synthesis = 0;
	double accuracy = 0.0;

	int L = 2;
	int P = 2;
	printf("  > B_l = %i - B_p = %i - J_min_l = %i - J_min_p = %i\n",B_l, B_p,J_min_l,J_min_p);

	for (sc=0; sc<NSCALE; sc++) {
		
		L *= 2;
		P *= 2;
	
		flag_allocate_flmn(&flmp, L, P);
		
		double *nodes = (double*)calloc(P, sizeof(double));
		double *weights = (double*)calloc(P, sizeof(double));
	 	flag_spherlaguerre_sampling(nodes, weights, R, P);

	 	printf("\n  > R = %2.2f - L = %i - P = %i \n", R, L, P);
	 	for (repeat=0; repeat<NREPEAT; repeat++){

	 		//printf("  -> Iteration : %i on %i\n",repeat+1,NREPEAT);

			flaglet_random_flmp(flmp, L, P, seed);

			flag_allocate_f(&f, L, P);
			flag_allocate_flmn(&flmp_rec, L, P);

			time_start = clock();
			flag_synthesis(f, flmp, nodes, P, L, P);
			time_end = clock();
			//printf("  - (FLAG synthesis   : %4.4f seconds)\n", (time_end - time_start) / (double)CLOCKS_PER_SEC);

			complex double *f_wav, *f_scal;
			flaglet_axisym_allocate_f_wav(&f_wav, &f_scal, B_l, B_p, L, P, J_min_l, J_min_p);

		    time_start = clock();
			flaglet_axisym_wav_analysis(f_wav, f_scal, f, R, B_l, B_p, L, P, J_min_l, J_min_p);
			time_end = clock();
			tottime_analysis += (time_end - time_start) / (double)CLOCKS_PER_SEC;
			//printf("  - Duration for wavelet analysis   : %4.4f seconds\n", (time_end - time_start) / (double)CLOCKS_PER_SEC);

			flag_allocate_f(&f_rec, L, P);

			time_start = clock();
			flaglet_axisym_wav_synthesis(f_rec, f_wav, f_scal, R, B_l, B_p, L, P, J_min_l, J_min_p);
			time_end = clock();
			tottime_synthesis += (time_end - time_start) / (double)CLOCKS_PER_SEC;
			//printf("  - Duration for wavelet synthesis   : %4.4f seconds\n", (time_end - time_start) / (double)CLOCKS_PER_SEC);
			
			time_start = clock();
			flag_analysis(flmp_rec, f_rec, R, L, P);
			time_end = clock();
			//printf("  - (FLAG analysis   : %4.4f seconds)\n", (time_end - time_start) / (double)CLOCKS_PER_SEC);

			accuracy += maxerr_cplx(flmp, flmp_rec, flag_flmn_size(L, P));
			//printf("  - Max error on reconstruction  : %6.5e\n\n", maxerr_cplx(flmp, flmp_rec, flag_flmn_size(L, P)));

			free(f);
			free(flmp_rec);
			free(f_rec);
			free(f_wav);
			free(f_scal);

		}

		tottime_analysis = tottime_analysis / (double)NREPEAT;
		tottime_synthesis = tottime_synthesis / (double)NREPEAT;
		accuracy = accuracy / (double)NREPEAT;
		
		printf("  - Average duration for wavelet analysis   : %4.4f seconds\n", tottime_analysis);
		printf("  - Average duration for wavelet synthesis  : %4.4f seconds\n", tottime_synthesis);
		printf("  - Average max error on reconstruction  : %6.5e\n", accuracy);

		free(flmp);
		free(nodes);
		free(weights);

	}

}

void flaglet_axisym_wav_multires_performance_test(double R, int NREPEAT, int NSCALE, int seed, int B_l, int B_p, int J_min_l, int J_min_p)
{
	complex double *f, *f_rec, *flmp, *flmp_rec;
	clock_t time_start, time_end;
	int sc, repeat;
	double tottime_analysis = 0, tottime_synthesis = 0;
	double accuracy = 0.0;

	int L = 2;
	int P = 2;
	printf("  > B_l = %i - B_p = %i - J_min_l = %i - J_min_p = %i\n",B_l, B_p,J_min_l,J_min_p);

	for (sc=0; sc<NSCALE; sc++) {
		
		L *= 2;
		P *= 2;
	
		flag_allocate_flmn(&flmp, L, P);
		
		double *nodes = (double*)calloc(P, sizeof(double));
		double *weights = (double*)calloc(P, sizeof(double));
	 	flag_spherlaguerre_sampling(nodes, weights, R, P);

	 	printf("\n  > R = %2.2f - L = %i - P = %i \n", R, L, P);
	 	for (repeat=0; repeat<NREPEAT; repeat++){

	 		//printf("  -> Iteration : %i on %i\n",repeat+1,NREPEAT);

			flaglet_random_flmp(flmp, L, P, seed);

			flag_allocate_f(&f, L, P);
			flag_allocate_flmn(&flmp_rec, L, P);

			time_start = clock();
			flag_synthesis(f, flmp, nodes, P, L, P);
			time_end = clock();
			//printf("  - (FLAG synthesis   : %4.4f seconds)\n", (time_end - time_start) / (double)CLOCKS_PER_SEC);

			complex double *f_wav, *f_scal;
			flaglet_axisym_allocate_f_wav_multires(&f_wav, &f_scal, B_l, B_p, L, P, J_min_l, J_min_p);

		    time_start = clock();
			flaglet_axisym_wav_analysis_multires(f_wav, f_scal, f, R, B_l, B_p, L, P, J_min_l, J_min_p);
			time_end = clock();
			tottime_analysis += (time_end - time_start) / (double)CLOCKS_PER_SEC;
			//printf("  - Duration for wavelet analysis   : %4.4f seconds\n", (time_end - time_start) / (double)CLOCKS_PER_SEC);

			flag_allocate_f(&f_rec, L, P);

			time_start = clock();
			flaglet_axisym_wav_synthesis_multires(f_rec, f_wav, f_scal, R, B_l, B_p, L, P, J_min_l, J_min_p);
			time_end = clock();
			tottime_synthesis += (time_end - time_start) / (double)CLOCKS_PER_SEC;
			//printf("  - Duration for wavelet synthesis   : %4.4f seconds\n", (time_end - time_start) / (double)CLOCKS_PER_SEC);
			
			time_start = clock();
			flag_analysis(flmp_rec, f_rec, R, L, P);
			time_end = clock();
			//printf("  - (FLAG analysis   : %4.4f seconds)\n", (time_end - time_start) / (double)CLOCKS_PER_SEC);

			accuracy += maxerr_cplx(flmp, flmp_rec, flag_flmn_size(L, P));
			//printf("  - Max error on reconstruction  : %6.5e\n\n", maxerr_cplx(flmp, flmp_rec, flag_flmn_size(L, P)));

			free(f);
			free(flmp_rec);
			free(f_rec);
			free(f_wav);
			free(f_scal);

		}

		tottime_analysis = tottime_analysis / (double)NREPEAT;
		tottime_synthesis = tottime_synthesis / (double)NREPEAT;
		accuracy = accuracy / (double)NREPEAT;
		
		printf("  - Average duration for wavelet analysis   : %4.4f seconds\n", tottime_analysis);
		printf("  - Average duration for wavelet synthesis  : %4.4f seconds\n", tottime_synthesis);
		printf("  - Average max error on reconstruction  : %6.5e\n", accuracy);

		free(flmp);
		free(nodes);
		free(weights);

	}

}

int main(int argc, char *argv[]) 
{
	const int NREPEAT = 4;
	const int NSCALE = 3;

	const double R = 1.0;
	const int L = 16;
	const int P = 16;
	const int B_l = 2;
	const int B_p = 2;
	const int J_min_l = 0;
	const int J_min_p = 0;
	const int seed = (int)(10000.0*(double)clock()/(double)CLOCKS_PER_SEC);

	printf("==========================================================\n");
	printf("PARAMETERS (seed = %i) \n", seed);
	printf(" L = %i   P = %i   Bl = %i   Bn = %i   Jminl = %i  Jminn = %i \n", L, P, B_l, B_p, J_min_l, J_min_p);
	printf("----------------------------------------------------------\n");
	printf("> Testing axisymmetric harmonic tilling...\n");
	flaglet_tilling_test(B_l, B_p, L, P, J_min_l, J_min_p);
	printf("==========================================================\n");
	printf("> Testing axisymmetric wavelets in harmonics space...\n");
	flaglet_axisym_wav_lm_test(B_l, B_p, L, P, J_min_l, J_min_p, seed);
	printf("----------------------------------------------------------\n");
	printf("> Testing multiresolution algorithm in harmonics space...\n");
	flaglet_axisym_wav_lm_multires_test(B_l, B_p, L, P, J_min_l, J_min_p, seed);
	printf("==========================================================\n");
	printf("> Testing axisymmetric wavelets in pixel space...\n");
	flaglet_axisym_wav_test(R, B_l, B_p, L, P, J_min_l, J_min_p, seed);
	printf("----------------------------------------------------------\n");
	printf("> Testing multiresolution algorithm...\n");
	flaglet_axisym_wav_multires_test(R, B_l, B_p, L, P, J_min_l, J_min_p, seed);
	printf("----------------------------------------------------------\n");
	printf("> Testing real axisymmetric wavelets in pixel space...\n");
	flaglet_axisym_wav_real_test(R, B_l, B_p, L, P, J_min_l, J_min_p, seed);
	printf("----------------------------------------------------------\n");
	printf("> Testing multiresolution algorithm...\n");
	flaglet_axisym_wav_multires_real_test(R, B_l, B_p, L, P, J_min_l, J_min_p, seed);
	printf("==========================================================\n");
	printf("> Full resolution wavelet transform : performance tests\n");
	flaglet_axisym_wav_performance_test(R, NREPEAT, NSCALE, seed, B_l, B_p, J_min_l, J_min_p);
	fflush(NULL);
	printf("----------------------------------------------------------\n");
	printf("> Multiresolution wavelet transform : performance tests\n");
	flaglet_axisym_wav_multires_performance_test(R, NREPEAT, NSCALE, seed, B_l, B_p, J_min_l, J_min_p);
	fflush(NULL);
	printf("==========================================================\n");
	
	return 0;		
}
