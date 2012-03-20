// B3LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "b3let.h"
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
    NTAB=32,NDIV=1+IMM1/NTAB;

  double AM=1./IM1,EPS=1.2e-7,RNMX=1.-EPS;
  int j,k;
  static int iv[32],iy,idum2 = 123456789; 
  // N.B. in C static variables are initialised to 0 by default.

  if (idum <= 0) {
    idum= (-idum>1 ? -idum : 1); // max(-idum,1);
    idum2=idum;
    for(j=NTAB+8;j>=1;j--) {
      k=idum/IQ1;
      idum=IA1*(idum-k*IQ1)-k*IR1;
      if (idum < 0) idum=idum+IM1;
      if (j < NTAB) iv[j-1]=idum;
    }
    iy=iv[0];
  }
  k=idum/IQ1;
  idum=IA1*(idum-k*IQ1)-k*IR1;
  if (idum < 0) idum=idum+IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2=idum2+IM2;
  j=1+iy/NDIV;
  iy=iv[j-1]-idum2;
  iv[j-1]=idum;
  if(iy < 1)iy=iy+IMM1;
  return (AM*iy < RNMX ? AM*iy : RNMX); // min(AM*iy,RNMX);
}

void b3let_random_flmn(complex double *flmn, int L, int N, int seed)
{
	int i;
	srand( time(NULL) );
	for (i=0; i<L*L*N; i++){
		flmn[i] = (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);
	}
}

void b3let_tilling_test(int B_l, int B_n, int L, int N, int J_min_l, int J_min_n)
{		
	double kl, kn;
	double *kappa_ln, *kappa0_ln;
	b3let_allocate_tilling(&kappa_ln, &kappa0_ln, B_l, B_n, L, N);
	
	b3let_tilling(kappa_ln, kappa0_ln, B_l, B_n, L, N, J_min_l, J_min_n);
	
	/*
	int jl, jn, l, n;
	int J_l = s2let_j_max(L, B_l);
	int J_n = s2let_j_max(N, B_n);
	printf("> KAPPA_0 :       ");
	for (n = 0; n < N; n++){
		for (l = 0; l < L; l++){
			printf(" %f ,", kappa0_ln[n * L + l]);
		}
		printf("\n");}
	printf("\n");
	for (jn = J_min_n; jn <= J_n; jn++){
		for (jl = J_min_l; jl <= J_l; jl++){
			for (n = 0; n < N; n++){
				for (l = 0; l < L; l++){
					printf(" %f ,",kappa_ln[jn*(J_l+1)*L*N  + jl*L*N + n*L + l]);
				}
				printf("\n");}
			printf("\n\n");}
		printf("\n\n\n");}
	*/

	double sum = b3let_check_identity(kappa_ln, kappa0_ln, B_l, B_n, L, N, J_min_l, J_min_n);
	printf("  - Identity residuals : %e\n", sum);

	free(kappa_ln);
	free(kappa0_ln);
}

void b3let_wav_lm_test(int B_l, int B_n, int L, int N, int J_min_l, int J_min_n, int seed)
{		
	clock_t time_start, time_end;

	double *wav_lmn, *scal_lmn;

	b3let_allocate_wav_lmn(&wav_lmn, &scal_lmn, B_l, B_n, L, N);

	time_start = clock();
	b3let_wav_lmn(wav_lmn, scal_lmn, B_l, B_n, L, N, J_min_l, J_min_n);
	time_end = clock();
	printf("  - Generate wavelets  : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	complex double *f_wav_lmn, *f_scal_lmn, *flmn, *flmn_rec;
	flmn = (complex double*)calloc(L * L * N, sizeof(complex double));
	flmn_rec = (complex double*)calloc(L * L * N, sizeof(complex double));

	b3let_random_flmn(flmn, L, N, seed);
	
	b3let_allocate_f_wav_lmn(&f_wav_lmn, &f_scal_lmn, B_l, B_n, L, N);

	time_start = clock();
	b3let_wav_analysis_lmn(f_wav_lmn, f_scal_lmn, flmn, wav_lmn, scal_lmn, B_l, B_n, L, N, J_min_l, J_min_n);
	time_end = clock();
	printf("  - Wavelet analysis   : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);


	int n, l, m;
	for (n = 0; n < N; n++){
		for (l = 0; l < L; l++){
	    	for (m = -l; m <= l ; m++){
		    	int ind = n*L*L + l*l + l + m;
				//printf("(l,m,n) = (%i,%i,%i) - f_scal_lmn = (%f,%f)\n",l,m,n,creal(f_scal_lmn[ind]),cimag(f_scal_lmn[ind]));
	   		}
	   	}
	}

	time_start = clock();
	b3let_wav_synthesis_lmn(flmn_rec, f_wav_lmn, f_scal_lmn, wav_lmn, scal_lmn, B_l, B_n, L, N, J_min_l, J_min_n);
	time_end = clock();
	printf("  - Wavelet synthesis  : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	for (n = 0; n < N; n++){
		for (l = 0; l < L; l++){
	    	for (m = -l; m <= l ; m++){
		    	int ind = n*L*L + l*l + l + m;
		    	if( creal(flmn[ind])-creal(flmn_rec[ind])>0.01 )
					printf("(l,m,n) = (%i,%i,%i) - flmn = (%f,%f) - rec = (%f,%f)\n",l,m,n,creal(flmn[ind]),cimag(flmn[ind]),creal(flmn_rec[ind]),cimag(flmn_rec[ind]));
	   		}
	   	}
	}

	printf("  - Maximum abs error  : %6.5e\n", 
		maxerr_cplx(flmn, flmn_rec, L*L*N));

	free(flmn);
	free(flmn_rec);
	free(f_wav_lmn);
	free(f_scal_lmn);
	free(wav_lmn);
	free(scal_lmn);

}


void b3let_wav_test(int B_l, int B_n, int L, int N, int J_min_l, int J_min_n, int seed)
{
	clock_t time_start, time_end;
	int l, m, n;
	int spin = 0;
	int verbosity = 0;
	int J_l = s2let_j_max(L, B_l);
	int J_n = s2let_j_max(N, B_n);

	complex double *f, *f_rec, *flmn, *flmn_rec;
	flag_allocate_flmn(&flmn, L, N);
	flag_allocate_flmn(&flmn_rec, L, N);
	flag_allocate_f(&f, L, N);
	flag_allocate_f(&f_rec, L, N);

	b3let_random_flmn(flmn, L, N, seed);

	flag_synthesis(f, flmn, L, N);

	complex double *f_wav, *f_scal;
	f_wav = (complex double*)calloc( (J_l+1) * L * (2*L-1) * (J_n+1) * N, sizeof(complex double));
	f_scal = (complex double*)calloc( L * (2*L-1) * N, sizeof(complex double));

	time_start = clock();
	b3let_wav_analysis(f_wav, f_scal, f, B_l, B_n, L, N, J_min_l, J_min_n);
	time_end = clock();
	printf("  - Wavelet analysis   : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	time_start = clock();
	b3let_wav_synthesis(f_rec, f_wav, f_scal, B_l, B_n, L, N, J_min_l, J_min_n);
	time_end = clock();
	printf("  - Wavelet synthesis  : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	/*printf("  - Maximum abs error  : %6.5e\n", 
		maxerr_cplx(f, f_rec, L*(2*L-1)*N));
	for (n = 0; n < L*(2*L-1)*N; n++){
		//printf("f[%i] = (%f,%f) - rec = (%f,%f)\n",n,creal(f[n]),cimag(f[n]),creal(f_rec[n]),cimag(f_rec[n]));
	}*/

	flag_analysis(flmn_rec, f_rec, L, N);

	printf("  - Maximum abs error  : %6.5e\n", 
		maxerr_cplx(flmn, flmn_rec, L*L*N));
	
	for (n = 0; n < N; n++){
		for (l = 0; l < L; l++){
			for (m = -l; m <= l ; m++){
		    	int ind = n*L*L+l*l+l+m;
				//printf("(l,m,n) = (%i,%i,%i) - flmn = (%f,%f) - rec = (%f,%f)\n",l,m,n,creal(flmn[ind]),cimag(flmn[ind]),creal(flmn_rec[ind]),cimag(flmn_rec[ind]));
	    	}
	    }
	}
	
	free(f);
	free(f_rec);
	free(f_wav);
	free(f_scal);
}

int main(int argc, char *argv[]) 
{

	const int L = 32;
	const int N = 32;
	const int B_l = 3;
	const int B_n = 2;
	const int J_min_l = 1;
	const int J_min_n = 3;
	const int seed = (int)(10000.0*(double)clock()/(double)CLOCKS_PER_SEC);

	printf("==========================================================\n");
	printf("PARAMETERS (seed = %i) \n", seed);
	printf(" L = %i   N = %i   Bl = %i   Bn = %i   Jminl = %i  Jminn = %i \n", L, N, B_l, B_n, J_min_l, J_min_n);
	printf("----------------------------------------------------------\n");
	printf("> Testing harmonic tilling...\n");
	b3let_tilling_test(B_l, B_n, L, N, J_min_l, J_min_n);
	printf("----------------------------------------------------------\n");
	printf("> Testing axisymmetric wavelets in harmonics space...\n");
	b3let_wav_lm_test(B_l, B_n, L, N, J_min_l, J_min_n, seed);
	printf("----------------------------------------------------------\n");
	printf("> Testing axisymmetric wavelets in pixel space...\n");
	b3let_wav_test(B_l, B_n, L, N, J_min_l, J_min_n, seed);
	printf("==========================================================\n");
	
	return 0;		
}
