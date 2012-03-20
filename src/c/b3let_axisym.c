// B3LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "b3let.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))

void b3let_allocate_wav_lmn(double **wav_lmn, double **scal_lmn, int B_l, int B_n, int L, int N)
{
	int J_l = s2let_j_max(L, B_l);
	int J_n = s2let_j_max(N, B_n);
	*wav_lmn = (double*)calloc( (J_l+1) * L * L * (J_n+1) * N, sizeof(double));
	*scal_lmn = (double*)calloc( L * L * N, sizeof(double));
}

void b3let_wav_lmn(double *wav_lmn, double *scal_lmn, int B_l, int B_n, int L, int N, int J_min_l, int J_min_n)
{
	int jl, jn, l, n, m, indjjlmn, indlmn;
	int J_l = s2let_j_max(L, B_l);
	int J_n = s2let_j_max(N, B_n);

	double wav0, scal0;

	double *kappa_ln, *kappa0_ln;
	b3let_allocate_tilling(&kappa_ln, &kappa0_ln, B_l, B_n, L, N);
	
	b3let_tilling(kappa_ln, kappa0_ln, B_l, B_n, L, N, J_min_l, J_min_n);

	for (jn = J_min_n; jn <= J_n; jn++){
		for (jl = J_min_l; jl <= J_l; jl++){
			for (n = 0; n < N; n++){
				for (l = 0; l < L; l++){
					wav0 = sqrt( (2 * l + 1) / (4.0 * PI) ) *
						kappa_ln[ jn*(J_l+1)*L*N  + jl*L*N + n*L + l ];
					for (m = -l; m <= l ; m++){
						indjjlmn = jjlmn2ind(jl,jn,l,m,n,J_l,J_n,L,N);
						wav_lmn[indjjlmn] =  wav0;
					}
				}
			}
		}
	}
	for (n = 0; n < N; n++){
		for (l = 0; l < L; l++){
			scal0 = sqrt( (2 * l + 1) / (4.0 * PI) ) *
				kappa0_ln[ n*L + l ];
			for (m = -l; m <= l ; m++){
				indlmn = lmn2ind(l,m,n,L);
				scal_lmn[indlmn] = scal0;
			}
		}
	}

	free(kappa_ln);
	free(kappa0_ln);
}

void b3let_allocate_f_wav_lmn(complex double **f_wav_lmn, complex double **f_scal_lmn, int B_l, int B_n, int L, int N)
{
	int J_l = s2let_j_max(L, B_l);
	int J_n = s2let_j_max(N, B_n);
	*f_wav_lmn = (complex double*)calloc( (J_l+1) * L * L * (J_n+1) * N, sizeof(complex double));
	*f_scal_lmn = (complex double*)calloc( L * L * N, sizeof(complex double));
}

void b3let_wav_analysis_lmn(complex double *f_wav_lmn, complex double *f_scal_lmn, const complex double *flmn, const double *wav_lmn, const double *scal_lmn, int B_l, int B_n, int L, int N, int J_min_l, int J_min_n)
{
	int jl, jn, l, m, n, indjjlmn, indlmn;
	int J_l = s2let_j_max(L, B_l);
	int J_n = s2let_j_max(N, B_n);
	double wav0, scal0;

	for (jn = J_min_n; jn <= J_n; jn++){
		for (jl = J_min_l; jl <= J_l; jl++){
			for (n = 0; n < N; n++){
				for (l = 0; l < L; l++){
					for (m = -l; m <= l ; m++){
						indjjlmn = jjlmn2ind(jl,jn,l,m,n,J_l,J_n,L,N);
						indlmn = lmn2ind(l,m,n,L);
						f_wav_lmn[indjjlmn] = flmn[indlmn] * sqrt((4.0*PI)/(2.0*l+1.0)) * wav_lmn[indjjlmn] ;
					}
				}
			}
		}
	}
	for (n = 0; n < N; n++){
		for (l = 0; l < L; l++){
			for (m = -l; m <= l ; m++){
				indlmn = lmn2ind(l,m,n,L);
				f_scal_lmn[indlmn] = flmn[indlmn] * sqrt((4.0*PI)/(2.0*l+1.0)) * scal_lmn[indlmn] ;
			}
		}
	}
}

void b3let_wav_synthesis_lmn(complex double *flmn, const complex double *f_wav_lmn, const complex double *f_scal_lmn, const double *wav_lmn, const double *scal_lmn, int B_l, int B_n, int L, int N, int J_min_l, int J_min_n)
{
	int jl, jn, l, m, n, indjjlmn, indlmn;
	int J_l = s2let_j_max(L, B_l);
	int J_n = s2let_j_max(N, B_n);
	double wav0, scal0;

	for (jn = J_min_n; jn <= J_n; jn++){
		for (jl = J_min_l; jl <= J_l; jl++){
			for (n = 0; n < N; n++){
				for (l = 0; l < L; l++){
					for (m = -l; m <= l ; m++){
						indjjlmn = jjlmn2ind(jl,jn,l,m,n,J_l,J_n,L,N);
						indlmn = lmn2ind(l,m,n,L);
						flmn[indlmn] += f_wav_lmn[indjjlmn] * sqrt((4.0*PI)/(2.0*l+1.0)) * wav_lmn[indjjlmn] ;
					}
				}
			}
		}
	}

	for (n = 0; n < N; n++){
		for (l = 0; l < L; l++){
			for (m = -l; m <= l ; m++){
				indlmn = lmn2ind(l,m,n,L);
				flmn[indlmn] += f_scal_lmn[indlmn] * sqrt((4.0*PI)/(2.0*l+1.0)) * scal_lmn[indlmn] ;
			}
		}
	}
	
}

void b3let_wav_analysis(complex double *f_wav, complex double *f_scal, const complex double *f, int B_l, int B_n, int L, int N, int J_min_l, int J_min_n)
{
	complex double *flmn;
	flag_allocate_flmn(&flmn, L, N);
	flag_analysis(flmn, f, L, N);

	complex double *wav_lmn, *scal_lmn, *f_wav_lmn, *f_scal_lmn;

	b3let_allocate_wav_lmn(&wav_lmn, &scal_lmn, B_l, B_n, L, N);
	b3let_wav_lmn(wav_lmn, scal_lmn, B_l, B_n, L, N, J_min_l, J_min_n);

	b3let_allocate_f_wav_lmn(&f_wav_lmn, &f_scal_lmn, B_l, B_n, L, N);
	b3let_wav_analysis_lmn(f_wav_lmn, f_scal_lmn, flmn, wav_lmn, scal_lmn, B_l, B_n, L, N, J_min_l, J_min_n);

	free(wav_lmn);
	free(scal_lmn);

	int offset_lmn, offset, jl, jn;
	int J_l = s2let_j_max(L, B_l);
	int J_n = s2let_j_max(N, B_n);
	flag_synthesis(f_scal, f_scal_lmn, L, N);
	for (jn = J_min_n; jn <= J_n; jn++){
		for (jl = J_min_l; jl <= J_l; jl++){
			offset_lmn = jn * (J_l + 1) * L * L * N   +  jl * L * L * N;
			offset = jn * (J_l + 1) * L * (2*L-1) * N   +  jl * L * (2*L-1) * N;
			flag_synthesis(f_wav + offset, f_wav_lmn + offset_lmn, L, N);
		}
	}

	free(f_wav_lmn);
	free(f_scal_lmn);
	free(flmn);
}

void b3let_wav_synthesis(complex double *f, const complex double *f_wav, const complex double *f_scal, int B_l, int B_n, int L, int N, int J_min_l, int J_min_n)
{
	complex double *f_wav_lmn, *f_scal_lmn;
	b3let_allocate_f_wav_lmn(&f_wav_lmn, &f_scal_lmn, B_l, B_n, L, N);

	int offset_lmn, offset, jl, jn;
	int J_l = s2let_j_max(L, B_l);
	int J_n = s2let_j_max(N, B_n);
	flag_analysis(f_scal_lmn, f_scal, L, N);
	for (jn = J_min_n; jn <= J_n; jn++){
		for (jl = J_min_l; jl <= J_l; jl++){
			offset_lmn = jn * (J_l + 1) * L * L * N   +  jl * L * L * N;
			offset = jn * (J_l + 1) * L * (2*L-1) * N   +  jl * L * (2*L-1) * N;
			flag_analysis(f_wav_lmn + offset_lmn, f_wav + offset, L, N);
		}
	}

	complex double *wav_lmn, *scal_lmn ;
	b3let_allocate_wav_lmn(&wav_lmn, &scal_lmn, B_l, B_n, L, N);
	b3let_wav_lmn(wav_lmn, scal_lmn, B_l, B_n, L, N, J_min_l, J_min_n);

	complex double *flmn;
	flag_allocate_flmn(&flmn, L, N);
	b3let_wav_synthesis_lmn(flmn, f_wav_lmn, f_scal_lmn, wav_lmn, scal_lmn, B_l, B_n, L, N, J_min_l, J_min_n);
	
	flag_synthesis(f, flmn, L, N);

	free(wav_lmn);
	free(scal_lmn);
	free(f_wav_lmn);
	free(f_scal_lmn);
	free(flmn);
}

int jjlmn2ind(int jl, int jn, int l, int m, int n, int J_l, int J_n, int L, int N)
{
	return jn * (J_l + 1) * L * L * N   +  jl * L * L * N  +  n * L * L  +  l * l  + l +  m;
}

int lmn2ind(int l, int m, int n, int L)
{
	return n * L * L  +  l * l  +  l  +  m ;
}


