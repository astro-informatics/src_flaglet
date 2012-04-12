// B3LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "b3let.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))

/*!
 * Allocates arrays for the kernels of the wavelets and the scaling functions (in FLAG space).
 *
 * \param[out]  wav_lmn Wavelet kernels.
 * \param[out]  scal_lmn Scaling function kernels.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_n Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  N Radial harmonic band-limit.
 * \retval none
 */
void b3let_axisym_allocate_wav_lmn(double **wav_lmn, double **scal_lmn, int B_l, int B_n, int L, int N)
{
	int J_l = s2let_j_max(L, B_l);
	int J_n = s2let_j_max(N, B_n);
	*wav_lmn = (double*)calloc( (J_l+1) * L * L * (J_n+1) * N, sizeof(double));
	*scal_lmn = (double*)calloc( L * L * N, sizeof(double));
}

/*!
 * Computes the kernels of the wavelets and the scaling functions (in FLAG space).
 *
 * \param[out]  wav_lmn Wavelet kernels.
 * \param[out]  scal_lmn Scaling function kernels.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_n Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  N Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_n First wavelet scale to be used in radial space.
 * \retval none
 */
void b3let_axisym_wav_lmn(double *wav_lmn, double *scal_lmn, int B_l, int B_n, int L, int N, int J_min_l, int J_min_n)
{
	int jl, jn, l, n, m, indjjlmn, indlmn;
	int J_l = s2let_j_max(L, B_l);
	int J_n = s2let_j_max(N, B_n);

	double wav0, scal0;

	double *kappa_ln, *kappa0_ln;
	b3let_axisym_allocate_tilling(&kappa_ln, &kappa0_ln, B_l, B_n, L, N);
	
	b3let_axisym_tilling(kappa_ln, kappa0_ln, B_l, B_n, L, N, J_min_l, J_min_n);

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

/*!
 * Allocates 3D Wavelet transform in FLAG space.
 *
 * \param[out]  f_wav_lmn FLAG transform of wavelets contributions of f.
 * \param[out]  f_scal_lmn FLAG transform of scaling function contribution of f.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_n Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  N Radial harmonic band-limit.
 * \retval none
 */
void b3let_axisym_allocate_f_wav_lmn(complex double **f_wav_lmn, complex double **f_scal_lmn, int B_l, int B_n, int L, int N)
{
	int J_l = s2let_j_max(L, B_l);
	int J_n = s2let_j_max(N, B_n);
	*f_wav_lmn = (complex double*)calloc( (J_l+1) * L * L * (J_n+1) * N, sizeof(complex double));
	*f_scal_lmn = (complex double*)calloc( L * L * N, sizeof(complex double));
}

/*!
 * Perform wavelet transform in FLAG space (from precomputed kernels, gives FLAG coefficients).
 * 3D spherical wavelets : analysis in FLAG-harmonic space.
 *
 * \param[out]  f_wav_lmn FLAG transform of wavelets contributions of f.
 * \param[out]  f_scal_lmn FLAG transform of scaling function contribution of f.
 * \param[in]  flmn FLAG transform of the input function.
 * \param[in]  wav_lmn Wavelet kernel in FLAG space.
 * \param[in]  scal_lmn Scaling function kernel in FLAG space.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_n Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  N Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_n First wavelet scale to be used in radial space.
 * \retval none
 */
void b3let_axisym_wav_analysis_lmn(complex double *f_wav_lmn, complex double *f_scal_lmn, const complex double *flmn, const double *wav_lmn, const double *scal_lmn, int B_l, int B_n, int L, int N, int J_min_l, int J_min_n)
{
	int jl, jn, l, m, n, indjjlmn, indlmn;
	int J_l = s2let_j_max(L, B_l);
	int J_n = s2let_j_max(N, B_n);

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

/*!
 * Perform wavelet transform in FLAG space (from precomputed kernels, gives FLAG coefficients).
 * 3D spherical wavelets : synthesis in FLAG-harmonic space.
 *
 * \param[out]  flmn FLAG transform of the input function.
 * \param[in]  f_wav_lmn FLAG transform of wavelets contributions of f.
 * \param[in]  f_scal_lmn FLAG transform of scaling function contribution of f.
 * \param[in]  wav_lmn Wavelet kernel in FLAG space.
 * \param[in]  scal_lmn Scaling function kernel in FLAG space.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_n Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  N Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_n First wavelet scale to be used in radial space.
 * \retval none
 */
void b3let_axisym_wav_synthesis_lmn(complex double *flmn, const complex double *f_wav_lmn, const complex double *f_scal_lmn, const double *wav_lmn, const double *scal_lmn, int B_l, int B_n, int L, int N, int J_min_l, int J_min_n)
{
	int jl, jn, l, m, n, indjjlmn, indlmn;
	int J_l = s2let_j_max(L, B_l);
	int J_n = s2let_j_max(N, B_n);

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

/*!
 * Perform wavelet transform in real space (from scratch, gives pixel space components).
 * Sampling scheme : MW sampling.
 * 3D spherical wavelets : analysis in real space.
 *
 * \param[out]  f_wav Wavelet transform (wavelet contribution in real space).
 * \param[out]  f_scal Wavelet transform (scaling contribution in real space).
 * \param[in]  f Input function (MW sampling).
 * \param[in]  wav_lmn Wavelet kernel in FLAG space.
 * \param[in]  scal_lmn Scaling function kernel in FLAG space.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_n Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  N Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_n First wavelet scale to be used in radial space.
 * \retval none
 */
void b3let_axisym_wav_analysis(complex double *f_wav, complex double *f_scal, const complex double *f, double R, int B_l, int B_n, int L, int N, int J_min_l, int J_min_n)
{
	complex double *flmn;
	flag_allocate_flmn(&flmn, L, N);
	flag_analysis(flmn, f, R, L, N);

	complex double *wav_lmn, *scal_lmn, *f_wav_lmn, *f_scal_lmn;

	b3let_axisym_allocate_wav_lmn(&wav_lmn, &scal_lmn, B_l, B_n, L, N);
	b3let_axisym_wav_lmn(wav_lmn, scal_lmn, B_l, B_n, L, N, J_min_l, J_min_n);

	b3let_axisym_allocate_f_wav_lmn(&f_wav_lmn, &f_scal_lmn, B_l, B_n, L, N);
	b3let_axisym_wav_analysis_lmn(f_wav_lmn, f_scal_lmn, flmn, wav_lmn, scal_lmn, B_l, B_n, L, N, J_min_l, J_min_n);

	free(wav_lmn);
	free(scal_lmn);

	double *nodes = (double*)calloc(N+1, sizeof(double));
	double *weights = (double*)calloc(N+1, sizeof(double));
	assert(nodes != NULL);
	assert(weights != NULL);
	flag_spherlaguerre_sampling(nodes, weights, R, N);

	int offset_lmn, offset, jl, jn;
	int J_l = s2let_j_max(L, B_l);
	int J_n = s2let_j_max(N, B_n);
	flag_synthesis(f_scal, f_scal_lmn, nodes, N+1, L, N);
	for (jn = J_min_n; jn <= J_n; jn++){
		for (jl = J_min_l; jl <= J_l; jl++){
			offset_lmn = jn * (J_l + 1) * L * L * N   +  jl * L * L * N;
			offset = jn * (J_l + 1) * L * (2*L-1) * (N+1)  +  jl * L * (2*L-1) * (N+1);
			flag_synthesis(f_wav + offset, f_wav_lmn + offset_lmn, nodes, N+1, L, N);
		}
	}

	free(nodes);
	free(weights);
	free(f_wav_lmn);
	free(f_scal_lmn);
	free(flmn);
}

/*!
 * Perform wavelet transform in real space (from scratch, gives pixel space components).
 * Sampling scheme : MW sampling.
 * 3D spherical wavelets : synthesis in real space.
 *
 * \param[out]  f Input function (MW sampling).
 * \param[in]  f_wav Wavelet transform (wavelet contribution in real space).
 * \param[in]  f_scal Wavelet transform (scaling contribution in real space).
 * \param[in]  wav_lmn Wavelet kernel in FLAG space.
 * \param[in]  scal_lmn Scaling function kernel in FLAG space.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_n Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  N Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_n First wavelet scale to be used in radial space.
 * \retval none
 */
void b3let_axisym_wav_synthesis(complex double *f, const complex double *f_wav, const complex double *f_scal, double R, int B_l, int B_n, int L, int N, int J_min_l, int J_min_n)
{
	complex double *f_wav_lmn, *f_scal_lmn;
	b3let_axisym_allocate_f_wav_lmn(&f_wav_lmn, &f_scal_lmn, B_l, B_n, L, N);

	int offset_lmn, offset, jl, jn;
	int J_l = s2let_j_max(L, B_l);
	int J_n = s2let_j_max(N, B_n);
	flag_analysis(f_scal_lmn, f_scal, R, L, N);
	for (jn = J_min_n; jn <= J_n; jn++){
		for (jl = J_min_l; jl <= J_l; jl++){
			offset_lmn = jn * (J_l + 1) * L * L * N   +  jl * L * L * N;
			offset = jn * (J_l + 1) * L * (2*L-1) * (N+1)   +  jl * L * (2*L-1) * (N+1);
			flag_analysis(f_wav_lmn + offset_lmn, f_wav + offset, R, L, N);
		}
	}

	complex double *wav_lmn, *scal_lmn ;
	b3let_axisym_allocate_wav_lmn(&wav_lmn, &scal_lmn, B_l, B_n, L, N);
	b3let_axisym_wav_lmn(wav_lmn, scal_lmn, B_l, B_n, L, N, J_min_l, J_min_n);

	complex double *flmn;
	flag_allocate_flmn(&flmn, L, N);
	b3let_axisym_wav_synthesis_lmn(flmn, f_wav_lmn, f_scal_lmn, wav_lmn, scal_lmn, B_l, B_n, L, N, J_min_l, J_min_n);
	
	double *nodes = (double*)calloc(N+1, sizeof(double));
	double *weights = (double*)calloc(N+1, sizeof(double));
	assert(nodes != NULL);
	assert(weights != NULL);
	flag_spherlaguerre_sampling(nodes, weights, R, N);

	flag_synthesis(f, flmn, nodes, N+1, L, N);

	free(nodes);
	free(weights);
	free(wav_lmn);
	free(scal_lmn);
	free(f_wav_lmn);
	free(f_scal_lmn);
	free(flmn);
}

/*!
 * Perform wavelet transform in real space (from scratch, gives pixel space components).
 * Real input function and real wavelet contributions.
 * Sampling scheme : MW sampling.
 * 3D spherical wavelets : analysis in real space.
 *
 * \param[out]  f_wav Wavelet transform (wavelet contribution in real space).
 * \param[out]  f_scal Wavelet transform (scaling contribution in real space).
 * \param[in]  f Input function (MW sampling).
 * \param[in]  wav_lmn Wavelet kernel in FLAG space.
 * \param[in]  scal_lmn Scaling function kernel in FLAG space.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_n Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  N Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_n First wavelet scale to be used in radial space.
 * \retval none
 */
void b3let_axisym_wav_analysis_real(double *f_wav, double *f_scal, const double *f, double R, int B_l, int B_n, int L, int N, int J_min_l, int J_min_n)
{
	complex double *flmn;
	flag_allocate_flmn(&flmn, L, N);
	flag_analysis_real(flmn, f, R, L, N);

	complex double *wav_lmn, *scal_lmn, *f_wav_lmn, *f_scal_lmn;

	b3let_axisym_allocate_wav_lmn(&wav_lmn, &scal_lmn, B_l, B_n, L, N);
	b3let_axisym_wav_lmn(wav_lmn, scal_lmn, B_l, B_n, L, N, J_min_l, J_min_n);

	b3let_axisym_allocate_f_wav_lmn(&f_wav_lmn, &f_scal_lmn, B_l, B_n, L, N);
	b3let_axisym_wav_analysis_lmn(f_wav_lmn, f_scal_lmn, flmn, wav_lmn, scal_lmn, B_l, B_n, L, N, J_min_l, J_min_n);

	free(wav_lmn);
	free(scal_lmn);

	double *nodes = (double*)calloc(N+1, sizeof(double));
	double *weights = (double*)calloc(N+1, sizeof(double));
	assert(nodes != NULL);
	assert(weights != NULL);
	flag_spherlaguerre_sampling(nodes, weights, R, N);

	int offset_lmn, offset, jl, jn;
	int J_l = s2let_j_max(L, B_l);
	int J_n = s2let_j_max(N, B_n);
	flag_synthesis_real(f_scal, f_scal_lmn, nodes, N+1, L, N);
	for (jn = J_min_n; jn <= J_n; jn++){
		for (jl = J_min_l; jl <= J_l; jl++){
			offset_lmn = jn * (J_l + 1) * L * L * N   +  jl * L * L * N;
			offset = jn * (J_l + 1) * L * (2*L-1) * (N+1)  +  jl * L * (2*L-1) * (N+1);
			flag_synthesis_real(f_wav + offset, f_wav_lmn + offset_lmn, nodes, N+1, L, N);
		}
	}

	free(nodes);
	free(weights);
	free(f_wav_lmn);
	free(f_scal_lmn);
	free(flmn);
}

/*!
 * Perform wavelet transform in real space (from scratch, gives pixel space components).
 * Real input function and real wavelet contributions.
 * Sampling scheme : MW sampling.
 * 3D spherical wavelets : synthesis in real space.
 *
 * \param[out]  f Input function (MW sampling).
 * \param[in]  f_wav Wavelet transform (wavelet contribution in real space).
 * \param[in]  f_scal Wavelet transform (scaling contribution in real space).
 * \param[in]  wav_lmn Wavelet kernel in FLAG space.
 * \param[in]  scal_lmn Scaling function kernel in FLAG space.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_n Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  N Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_n First wavelet scale to be used in radial space.
 * \retval none
 */
void b3let_axisym_wav_synthesis_real(double *f, const double *f_wav, const double *f_scal, double R, int B_l, int B_n, int L, int N, int J_min_l, int J_min_n)
{
	complex double *f_wav_lmn, *f_scal_lmn;
	b3let_axisym_allocate_f_wav_lmn(&f_wav_lmn, &f_scal_lmn, B_l, B_n, L, N);

	int offset_lmn, offset, jl, jn;
	int J_l = s2let_j_max(L, B_l);
	int J_n = s2let_j_max(N, B_n);
	flag_analysis_real(f_scal_lmn, f_scal, R, L, N);
	for (jn = J_min_n; jn <= J_n; jn++){
		for (jl = J_min_l; jl <= J_l; jl++){
			offset_lmn = jn * (J_l + 1) * L * L * N   +  jl * L * L * N;
			offset = jn * (J_l + 1) * L * (2*L-1) * (N+1)   +  jl * L * (2*L-1) * (N+1);
			flag_analysis_real(f_wav_lmn + offset_lmn, f_wav + offset, R, L, N);
		}
	}

	complex double *wav_lmn, *scal_lmn ;
	b3let_axisym_allocate_wav_lmn(&wav_lmn, &scal_lmn, B_l, B_n, L, N);
	b3let_axisym_wav_lmn(wav_lmn, scal_lmn, B_l, B_n, L, N, J_min_l, J_min_n);

	complex double *flmn;
	flag_allocate_flmn(&flmn, L, N);
	b3let_axisym_wav_synthesis_lmn(flmn, f_wav_lmn, f_scal_lmn, wav_lmn, scal_lmn, B_l, B_n, L, N, J_min_l, J_min_n);
	
	double *nodes = (double*)calloc(N+1, sizeof(double));
	double *weights = (double*)calloc(N+1, sizeof(double));
	assert(nodes != NULL);
	assert(weights != NULL);
	flag_spherlaguerre_sampling(nodes, weights, R, N);

	flag_synthesis_real(f, flmn, nodes, N+1, L, N);

	free(nodes);
	free(weights);
	free(wav_lmn);
	free(scal_lmn);
	free(f_wav_lmn);
	free(f_scal_lmn);
	free(flmn);
}

/*!
 * Indice corresponding to a quintuplet (jl, jn, l, m, n) in the wavelets.
 *
 * \param[in]  jl Angular scale indice.
 * \param[in]  jn Radial scale indice.
 * \param[in]  l Multipole indice.
 * \param[in]  m Order indice.
 * \param[in]  n Laguerre order indice.
 * \param[in]  J_l Maximum scale for angular harmonic space.
 * \param[in]  J_n Maximum scale for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  N Radial harmonic band-limit.
 * \retval Indice
 */
int jjlmn2ind(int jl, int jn, int l, int m, int n, int J_l, int J_n, int L, int N)
{
	return jn * (J_l + 1) * L * L * N   +  jl * L * L * N  +  n * L * L  +  l * l  + l +  m;
}

/*!
 * Indice corresponding to a triplet (l, m, n) in the FLAG basis.
 *
 * \param[in]  l Multipole indice.
 * \param[in]  m Order indice.
 * \param[in]  n Laguerre order indice.
 * \param[in]  L Angular harmonic band-limit.
 * \retval Indice
 */
int lmn2ind(int l, int m, int n, int L)
{
	return n * L * L  +  l * l  +  l  +  m ;
}


