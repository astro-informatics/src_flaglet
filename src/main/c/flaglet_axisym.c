// FLAGLET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "flaglet.h"
#include <stdlib.h>
#include <math.h>
#include <complex.h> 
#include <flag.h>
#include <s2let.h>
#include <assert.h>

/*!
 * Allocates arrays for the kernels of the wavelets and the scaling functions (in FLAG space).
 *
 * \param[out]  wav_lmp Wavelet kernels.
 * \param[out]  scal_lmp Scaling function kernels.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_p Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  P Radial harmonic band-limit.
 * \retval none
 */
void flaglet_axisym_allocate_wav_lmp(double **wav_lmp, double **scal_lmp, int B_l, int B_p, int L, int P)
{
	int J_l = s2let_j_max(L, B_l);
	int J_p = s2let_j_max(P, B_p);
	*wav_lmp = (double*)calloc( (J_l+1) * L * L * (J_p+1) * P, sizeof(double));
	*scal_lmp = (double*)calloc( L * L * P, sizeof(double));
}

/*!
 * Computes the kernels of the wavelets and the scaling functions (in FLAG space).
 *
 * \param[out]  wav_lmp Wavelet kernels.
 * \param[out]  scal_lmp Scaling function kernels.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_p Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  P Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_p First wavelet scale to be used in radial space.
 * \retval none
 */
void flaglet_axisym_wav_lmp(double *wav_lmp, double *scal_lmp, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p)
{
	int jl, jp, l, n, m, indjjlmp, indlmp;
	int J_l = s2let_j_max(L, B_l);
	int J_p = s2let_j_max(P, B_p);

	double wav0, scal0;

	double *kappa_ln, *kappa0_ln;
	flaglet_tiling_axisym_allocate(&kappa_ln, &kappa0_ln, B_l, B_p, L, P);
	
	flaglet_tiling_axisym(kappa_ln, kappa0_ln, B_l, B_p, L, P, J_min_l, J_min_p);

	for (jp = J_min_p; jp <= J_p; jp++){
		for (jl = J_min_l; jl <= J_l; jl++){
			for (n = 0; n < P; n++){
				for (l = 0; l < L; l++){
					wav0 = sqrt( (2 * l + 1) / (4.0 * PI) ) *
						kappa_ln[ jp*(J_l+1)*L*P  + jl*L*P + n*L + l ];
					for (m = -l; m <= l ; m++){
						indjjlmp = jjlmp2ind(jl,jp,l,m,n,J_l,J_p,L,P);
						wav_lmp[indjjlmp] =  wav0;
					}
				}
			}
		}
	}
	for (n = 0; n < P; n++){
		for (l = 0; l < L; l++){
			scal0 = sqrt( (2 * l + 1) / (4.0 * PI) ) *
				kappa0_ln[ n*L + l ];
			for (m = -l; m <= l ; m++){
				indlmp = lmp2ind(l,m,n,L);
				scal_lmp[indlmp] = scal0;
			}
		}
	}

	free(kappa_ln);
	free(kappa0_ln);
}

/*!
 * Allocates 3D Wavelet transform in FLAG space.
 *
 * \param[out]  f_wav_lmp FLAG transform of wavelets contributions of f.
 * \param[out]  f_scal_lmp FLAG transform of scaling function contribution of f.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_p Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  P Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_p First wavelet scale to be used in radial space.
 * \retval none
 */
void flaglet_axisym_allocate_f_wav_lmp(complex double **f_wav_lmp, complex double **f_scal_lmp, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p)
{
	int J_l = s2let_j_max(L, B_l);
	int J_p = s2let_j_max(P, B_p);
	*f_wav_lmp = (complex double*)calloc( (J_l+1-J_min_l) * L * L * (J_p+1-J_min_p) * P, sizeof(complex double));
	*f_scal_lmp = (complex double*)calloc( L * L * P, sizeof(complex double));
}

/*!
 * Allocates 3D Wavelet transform in real space.
 *
 * \param[out]  f_wav wavelets contributions of f.
 * \param[out]  f_scal scaling function contribution of f.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_p Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  P Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_p First wavelet scale to be used in radial space.
 * \retval none
 */
void flaglet_axisym_allocate_f_wav(complex double **f_wav, complex double **f_scal, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p)
{
	int J_l = s2let_j_max(L, B_l);
	int J_p = s2let_j_max(P, B_p);
	*f_wav = (complex double*)calloc( (J_l+1-J_min_l) * L * (2*L-1) * (J_p+1-J_min_p) * P, sizeof(complex double));
	*f_scal = (complex double*)calloc( L * (2*L-1) * P, sizeof(complex double));
}

/*!
 * Allocates real 3D Wavelet transform in real space.
 *
 * \param[out]  f_wav wavelets contributions of f.
 * \param[out]  f_scal scaling function contribution of f.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_p Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  P Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_p First wavelet scale to be used in radial space.
 * \retval none
 */
void flaglet_axisym_allocate_f_wav_real(double **f_wav, double **f_scal, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p)
{
	int J_l = s2let_j_max(L, B_l);
	int J_p = s2let_j_max(P, B_p);
	*f_wav = (double*)calloc( (J_l+1-J_min_l) * L * (2*L-1) * (J_p+1-J_min_p) * P, sizeof(double));
	*f_scal = (double*)calloc( L * (2*L-1) * P, sizeof(double));
}

/*!
 * Allocates multiresolution 3D Wavelet transform in FLAG space.
 *
 * \param[out]  f_wav_lmp FLAG transform of wavelets contributions of f.
 * \param[out]  f_scal_lmp FLAG transform of scaling function contribution of f.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_p Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  P Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_p First wavelet scale to be used in radial space.
 * \retval none
 */
void flaglet_axisym_allocate_f_wav_multires_lmp(complex double **f_wav_lmp, complex double **f_scal_lmp, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p)
{
	int J_l = s2let_j_max(L, B_l);
	int J_p = s2let_j_max(P, B_p);
	int jp, jl, bandlimit_p, bandlimit_l, total = 0;
	for (jp = J_min_p; jp <= J_p; jp++){
		bandlimit_p = MIN(s2let_bandlimit(B_p, jp), P);
		for (jl = J_min_l; jl <= J_l; jl++){
			bandlimit_l = MIN(s2let_bandlimit(B_l, jl), L);
			total += bandlimit_l * bandlimit_l * bandlimit_p;
		}
	}
	*f_wav_lmp = (complex double*)calloc( total, sizeof(complex double));
	*f_scal_lmp = (complex double*)calloc( L * L * P, sizeof(complex double));
}

/*!
 * Allocates multiresolution 3D Wavelet transform in real space.
 *
 * \param[out]  f_wav wavelets contributions of f.
 * \param[out]  f_scal scaling function contribution of f.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_p Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  P Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_p First wavelet scale to be used in radial space.
 * \retval none
 */
void flaglet_axisym_allocate_f_wav_multires(complex double **f_wav, complex double **f_scal, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p)
{
	int J_l = s2let_j_max(L, B_l);
	int J_p = s2let_j_max(P, B_p);
	int jp, jl, bandlimit_p, bandlimit_l, total = 0;
	for (jp = J_min_p; jp <= J_p; jp++){
		bandlimit_p = MIN(s2let_bandlimit(B_p, jp), P);
		for (jl = J_min_l; jl <= J_l; jl++){
			bandlimit_l = MIN(s2let_bandlimit(B_l, jl), L);
			total += bandlimit_l * (2 * bandlimit_l - 1) * bandlimit_p;
		}
	}
	*f_wav = (complex double*)calloc( total, sizeof(complex double));
	*f_scal = (complex double*)calloc( L * (2*L-1) * P, sizeof(complex double));
}

/*!
 * Allocates multiresolution real 3D Wavelet transform in real space.
 *
 * \param[out]  f_wav wavelets contributions of f.
 * \param[out]  f_scal scaling function contribution of f.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_p Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  P Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_p First wavelet scale to be used in radial space.
 * \retval none
 */
void flaglet_axisym_allocate_f_wav_multires_real(double **f_wav, double **f_scal, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p)
{
	int J_l = s2let_j_max(L, B_l);
	int J_p = s2let_j_max(P, B_p);
	int jp, jl, bandlimit_p, bandlimit_l, total = 0;
	for (jp = J_min_p; jp <= J_p; jp++){
		bandlimit_p = MIN(s2let_bandlimit(B_p, jp), P);
		for (jl = J_min_l; jl <= J_l; jl++){
			bandlimit_l = MIN(s2let_bandlimit(B_l, jl), L);
			total += bandlimit_l * (2 * bandlimit_l - 1) * bandlimit_p;
		}
	}
	*f_wav = (double*)calloc( total, sizeof(double));
	*f_scal = (double*)calloc( L * (2*L-1) * P, sizeof(double));
}

/*!
 * Perform multiresolution wavelet transform in FLAG space (from precomputed kernels, gives FLAG coefficients).
 * 3D spherical wavelets : analysis in FLAG-harmonic space.
 *
 * \param[out]  f_wav_lmp FLAG transform of wavelets contributions of f.
 * \param[out]  f_scal_lmp FLAG transform of scaling function contribution of f.
 * \param[in]  flmp FLAG transform of the input function.
 * \param[in]  wav_lmp Wavelet kernel in FLAG space.
 * \param[in]  scal_lmp Scaling function kernel in FLAG space.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_p Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  P Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_p First wavelet scale to be used in radial space.
 * \retval none
 */
void flaglet_axisym_wav_analysis_multires_lmp(complex double *f_wav_lmp, complex double *f_scal_lmp, const complex double *flmp, const double *wav_lmp, const double *scal_lmp, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p)
{
	int bandlimit_p, bandlimit_l, offset, jl, jp, l, m, n, indjjlmp, indlmp;
	int J_l = s2let_j_max(L, B_l);
	int J_p = s2let_j_max(P, B_p);

	offset = 0;
	for (jp = J_min_p; jp <= J_p; jp++){
		bandlimit_p = MIN(s2let_bandlimit(B_p, jp), P);
		for (jl = J_min_l; jl <= J_l; jl++){
			bandlimit_l = MIN(s2let_bandlimit(B_l, jl), L);
			for (n = 0; n < bandlimit_p; n++){
				for (l = 0; l < bandlimit_l; l++){
					for (m = -l; m <= l ; m++){
						indjjlmp = jjlmp2ind(jl,jp,l,m,n,J_l,J_p,L,P);
						indlmp = lmp2ind(l,m,n,L);
						f_wav_lmp[offset + lmp2ind(l,m,n,bandlimit_l)] = flmp[indlmp] * sqrt((4.0*PI)/(2.0*l+1.0)) * wav_lmp[indjjlmp] ;
					}
				}
			}
			offset += bandlimit_l * bandlimit_l * bandlimit_p;
		}
	}

	for (n = 0; n < P; n++){
		for (l = 0; l < L; l++){
			for (m = -l; m <= l ; m++){
				indlmp = lmp2ind(l,m,n,L);
				f_scal_lmp[indlmp] = flmp[indlmp] * sqrt((4.0*PI)/(2.0*l+1.0)) * scal_lmp[indlmp] ;
			}
		}
	}
}

/*!
 * Perform multiresolution wavelet transform in FLAG space (from precomputed kernels, gives FLAG coefficients).
 * 3D spherical wavelets : synthesis in FLAG-harmonic space.
 *
 * \param[out]  flmp FLAG transform of the input function.
 * \param[in]  f_wav_lmp FLAG transform of wavelets contributions of f.
 * \param[in]  f_scal_lmp FLAG transform of scaling function contribution of f.
 * \param[in]  wav_lmp Wavelet kernel in FLAG space.
 * \param[in]  scal_lmp Scaling function kernel in FLAG space.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_p Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  P Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_p First wavelet scale to be used in radial space.
 * \retval none
 */
void flaglet_axisym_wav_synthesis_multires_lmp(complex double *flmp, const complex double *f_wav_lmp, const complex double *f_scal_lmp, const double *wav_lmp, const double *scal_lmp, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p)
{
	int bandlimit_p, bandlimit_l, offset, jl, jp, l, m, n, indjjlmp, indlmp;
	int J_l = s2let_j_max(L, B_l);
	int J_p = s2let_j_max(P, B_p);

	offset = 0;
	for (jp = J_min_p; jp <= J_p; jp++){
		bandlimit_p = MIN(s2let_bandlimit(B_p, jp), P);
		for (jl = J_min_l; jl <= J_l; jl++){
			bandlimit_l = MIN(s2let_bandlimit(B_l, jl), L);
			for (n = 0; n < bandlimit_p; n++){
				for (l = 0; l < bandlimit_l; l++){
					for (m = -l; m <= l ; m++){
						indjjlmp = jjlmp2ind(jl,jp,l,m,n,J_l,J_p,L,P);
						indlmp = lmp2ind(l,m,n,L);
						flmp[indlmp] += f_wav_lmp[offset + lmp2ind(l,m,n,bandlimit_l)] * sqrt((4.0*PI)/(2.0*l+1.0)) * wav_lmp[indjjlmp] ;
					}
				}
			}
			offset += bandlimit_l * bandlimit_l * bandlimit_p;
		}
	}

	bandlimit_p = MIN(s2let_bandlimit(B_p, J_min_p-1), P);
	bandlimit_l = MIN(s2let_bandlimit(B_l, J_min_l-1), L);
	for (n = 0; n < P; n++){
		for (l = 0; l < L; l++){
			for (m = -l; m <= l ; m++){
				indlmp = lmp2ind(l,m,n,L);
				flmp[indlmp] += f_scal_lmp[indlmp] * sqrt((4.0*PI)/(2.0*l+1.0)) * scal_lmp[indlmp] ;
			}
		}
	}
	
}

/*!
 * Perform wavelet transform in FLAG space (from precomputed kernels, gives FLAG coefficients).
 * 3D spherical wavelets : analysis in FLAG-harmonic space.
 *
 * \param[out]  f_wav_lmp FLAG transform of wavelets contributions of f.
 * \param[out]  f_scal_lmp FLAG transform of scaling function contribution of f.
 * \param[in]  flmp FLAG transform of the input function.
 * \param[in]  wav_lmp Wavelet kernel in FLAG space.
 * \param[in]  scal_lmp Scaling function kernel in FLAG space.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_p Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  P Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_p First wavelet scale to be used in radial space.
 * \retval none
 */
void flaglet_axisym_wav_analysis_lmp(complex double *f_wav_lmp, complex double *f_scal_lmp, const complex double *flmp, const double *wav_lmp, const double *scal_lmp, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p)
{
	int offset, jl, jp, l, m, n, indjjlmp, indlmp;
	int J_l = s2let_j_max(L, B_l);
	int J_p = s2let_j_max(P, B_p);

	offset = 0;
	for (jp = J_min_p; jp <= J_p; jp++){
		for (jl = J_min_l; jl <= J_l; jl++){
			for (n = 0; n < P; n++){
				for (l = 0; l < L; l++){
					for (m = -l; m <= l ; m++){
						indjjlmp = jjlmp2ind(jl,jp,l,m,n,J_l,J_p,L,P);
						indlmp = lmp2ind(l,m,n,L);
						f_wav_lmp[offset + indlmp] = flmp[indlmp] * sqrt((4.0*PI)/(2.0*l+1.0)) * wav_lmp[indjjlmp] ;
					}
				}
			}
			offset += L * L * P;
		}
	}
	for (n = 0; n < P; n++){
		for (l = 0; l < L; l++){
			for (m = -l; m <= l ; m++){
				indlmp = lmp2ind(l,m,n,L);
				f_scal_lmp[indlmp] = flmp[indlmp] * sqrt((4.0*PI)/(2.0*l+1.0)) * scal_lmp[indlmp] ;
			}
		}
	}
}

/*!
 * Perform wavelet transform in FLAG space (from precomputed kernels, gives FLAG coefficients).
 * 3D spherical wavelets : synthesis in FLAG-harmonic space.
 *
 * \param[out]  flmp FLAG transform of the input function.
 * \param[in]  f_wav_lmp FLAG transform of wavelets contributions of f.
 * \param[in]  f_scal_lmp FLAG transform of scaling function contribution of f.
 * \param[in]  wav_lmp Wavelet kernel in FLAG space.
 * \param[in]  scal_lmp Scaling function kernel in FLAG space.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_p Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  P Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_p First wavelet scale to be used in radial space.
 * \retval none
 */
void flaglet_axisym_wav_synthesis_lmp(complex double *flmp, const complex double *f_wav_lmp, const complex double *f_scal_lmp, const double *wav_lmp, const double *scal_lmp, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p)
{
	int offset, jl, jp, l, m, n, indjjlmp, indlmp;
	int J_l = s2let_j_max(L, B_l);
	int J_p = s2let_j_max(P, B_p);

	offset = 0;
	for (jp = J_min_p; jp <= J_p; jp++){
		for (jl = J_min_l; jl <= J_l; jl++){
			for (n = 0; n < P; n++){
				for (l = 0; l < L; l++){
					for (m = -l; m <= l ; m++){
						indjjlmp = jjlmp2ind(jl,jp,l,m,n,J_l,J_p,L,P);
						indlmp = lmp2ind(l,m,n,L);
						flmp[indlmp] += f_wav_lmp[offset + indlmp] * sqrt((4.0*PI)/(2.0*l+1.0)) * wav_lmp[indjjlmp] ;
					}
				}
			}
			offset += L * L * P;
		}
	}

	for (n = 0; n < P; n++){
		for (l = 0; l < L; l++){
			for (m = -l; m <= l ; m++){
				indlmp = lmp2ind(l,m,n,L);
				flmp[indlmp] += f_scal_lmp[indlmp] * sqrt((4.0*PI)/(2.0*l+1.0)) * scal_lmp[indlmp] ;
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
 * \param[in]  wav_lmp Wavelet kernel in FLAG space.
 * \param[in]  scal_lmp Scaling function kernel in FLAG space.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_p Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  P Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_p First wavelet scale to be used in radial space.
 * \retval none
 */
void flaglet_axisym_wav_analysis_multires(complex double *f_wav, complex double *f_scal, const complex double *f, double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p)
{
	complex double *flmp;
	flag_core_allocate_flmn(&flmp, L, P);
	flag_core_analysis(flmp, f, R, L, P);

	double *wav_lmp, *scal_lmp;
	flaglet_axisym_allocate_wav_lmp(&wav_lmp, &scal_lmp, B_l, B_p, L, P);
	flaglet_axisym_wav_lmp(wav_lmp, scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);

	complex double *f_wav_lmp, *f_scal_lmp;
	flaglet_axisym_allocate_f_wav_multires_lmp(&f_wav_lmp, &f_scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);
	flaglet_axisym_wav_analysis_multires_lmp(f_wav_lmp, f_scal_lmp, flmp, wav_lmp, scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);

	double *nodes = (double*)calloc(P, sizeof(double));
	double *weights = (double*)calloc(P, sizeof(double));
	assert(nodes != NULL);
	assert(weights != NULL);
	flag_spherlaguerre_sampling(nodes, weights, R, P);

	int bandlimit_p, bandlimit_l, offset_lmp, offset, jl, jp;
	int J_l = s2let_j_max(L, B_l);
	int J_p = s2let_j_max(P, B_p);

	flag_core_synthesis(f_scal, f_scal_lmp, nodes, P, L, P);

	offset_lmp = 0;
	offset = 0;
	for (jp = J_min_p; jp <= J_p; jp++){
		bandlimit_p = MIN(s2let_bandlimit(B_p, jp), P);
		flag_spherlaguerre_sampling(nodes, weights, R, bandlimit_p);
		for (jl = J_min_l; jl <= J_l; jl++){
			bandlimit_l = MIN(s2let_bandlimit(B_l, jl), L);			
			//offset_lmp = jp * (J_l + 1) * L * L * P   +  jl * L * L * P;
			//offset = jp * (J_l + 1) * L * (2*L-1) * (P)  +  jl * L * (2*L-1) * (P);
			flag_core_synthesis(f_wav + offset, f_wav_lmp + offset_lmp, nodes, bandlimit_p, bandlimit_l, bandlimit_p);
			offset_lmp += bandlimit_l * bandlimit_l * bandlimit_p ;
			offset += bandlimit_l * (2 * bandlimit_l - 1) * bandlimit_p;
		}
	}

	free(nodes);
	free(weights);
	free(f_wav_lmp);
	free(f_scal_lmp);
	free(flmp);
}

/*!
 * Perform wavelet transform in real space (from scratch, gives pixel space components).
 * Sampling scheme : MW sampling.
 * 3D spherical wavelets : synthesis in real space.
 *
 * \param[out]  f Input function (MW sampling).
 * \param[in]  f_wav Wavelet transform (wavelet contribution in real space).
 * \param[in]  f_scal Wavelet transform (scaling contribution in real space).
 * \param[in]  wav_lmp Wavelet kernel in FLAG space.
 * \param[in]  scal_lmp Scaling function kernel in FLAG space.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_p Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  P Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_p First wavelet scale to be used in radial space.
 * \retval none
 */
void flaglet_axisym_wav_synthesis_multires(complex double *f, const complex double *f_wav, const complex double *f_scal, double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p)
{
	complex double *f_wav_lmp, *f_scal_lmp;
	flaglet_axisym_allocate_f_wav_multires_lmp(&f_wav_lmp, &f_scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);

	int bandlimit_p, bandlimit_l, offset_lmp, offset, jl, jp;
	int J_l = s2let_j_max(L, B_l);
	int J_p = s2let_j_max(P, B_p);

	flag_core_analysis(f_scal_lmp, f_scal, R, L, P);

	offset_lmp = 0;
	offset = 0;
	for (jp = J_min_p; jp <= J_p; jp++){
		bandlimit_p = MIN(s2let_bandlimit(B_p, jp), P);
		for (jl = J_min_l; jl <= J_l; jl++){
			bandlimit_l = MIN(s2let_bandlimit(B_l, jl), L);
			//offset_lmp = jp * (J_l + 1) * L * L * P   +  jl * L * L * P;
			//offset = jp * (J_l + 1) * L * (2*L-1) * (P)   +  jl * L * (2*L-1) * (P);
			flag_core_analysis(f_wav_lmp + offset_lmp, f_wav + offset, R, bandlimit_l, bandlimit_p);
			offset_lmp += bandlimit_l * bandlimit_l * bandlimit_p ;
			offset += bandlimit_l * (2 * bandlimit_l - 1) * bandlimit_p;
		}
	}

	double *wav_lmp, *scal_lmp ;
	flaglet_axisym_allocate_wav_lmp(&wav_lmp, &scal_lmp, B_l, B_p, L, P);
	flaglet_axisym_wav_lmp(wav_lmp, scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);

	complex double *flmp;
	flag_core_allocate_flmn(&flmp, L, P);
	flaglet_axisym_wav_synthesis_multires_lmp(flmp, f_wav_lmp, f_scal_lmp, wav_lmp, scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);


	double *nodes = (double*)calloc(P, sizeof(double));
	double *weights = (double*)calloc(P, sizeof(double));
	assert(nodes != NULL);
	assert(weights != NULL);
	flag_spherlaguerre_sampling(nodes, weights, R, P);

	flag_core_synthesis(f, flmp, nodes, P, L, P);

	free(nodes);
	free(weights);
	free(wav_lmp);
	free(scal_lmp);
	free(f_wav_lmp);
	free(f_scal_lmp);
	free(flmp);
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
 * \param[in]  wav_lmp Wavelet kernel in FLAG space.
 * \param[in]  scal_lmp Scaling function kernel in FLAG space.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_p Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  P Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_p First wavelet scale to be used in radial space.
 * \retval none
 */
void flaglet_axisym_wav_analysis_multires_real(double *f_wav, double *f_scal, const double *f, double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p)
{
	complex double *flmp;
	flag_core_allocate_flmn(&flmp, L, P);
	flag_core_analysis_real(flmp, f, R, L, P);

	double *wav_lmp, *scal_lmp;
	complex double *f_wav_lmp, *f_scal_lmp;

	flaglet_axisym_allocate_wav_lmp(&wav_lmp, &scal_lmp, B_l, B_p, L, P);
	flaglet_axisym_wav_lmp(wav_lmp, scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);

	flaglet_axisym_allocate_f_wav_multires_lmp(&f_wav_lmp, &f_scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);
	flaglet_axisym_wav_analysis_multires_lmp(f_wav_lmp, f_scal_lmp, flmp, wav_lmp, scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);

	free(wav_lmp);
	free(scal_lmp);

	double *nodes = (double*)calloc(P, sizeof(double));
	double *weights = (double*)calloc(P, sizeof(double));
	assert(nodes != NULL);
	assert(weights != NULL);
	flag_spherlaguerre_sampling(nodes, weights, R, P);

	int bandlimit_p, bandlimit_l, offset_lmp, offset, jl, jp;
	int J_l = s2let_j_max(L, B_l);
	int J_p = s2let_j_max(P, B_p);
	flag_core_synthesis_real(f_scal, f_scal_lmp, nodes, P, L, P);

	offset_lmp = 0;
	offset = 0;
	for (jp = J_min_p; jp <= J_p; jp++){
		bandlimit_p = MIN(s2let_bandlimit(B_p, jp), P);
		flag_spherlaguerre_sampling(nodes, weights, R, bandlimit_p);
		for (jl = J_min_l; jl <= J_l; jl++){
			bandlimit_l = MIN(s2let_bandlimit(B_l, jl), L);
			//offset_lmp = jp * (J_l + 1) * L * L * P   +  jl * L * L * P;
			//offset = jp * (J_l + 1) * L * (2*L-1) * (P)  +  jl * L * (2*L-1) * (P);
			flag_core_synthesis_real(f_wav + offset, f_wav_lmp + offset_lmp, nodes, bandlimit_p, bandlimit_l, bandlimit_p);
			offset_lmp += bandlimit_l * bandlimit_l * bandlimit_p ;
			offset += bandlimit_l * (2 * bandlimit_l - 1) * bandlimit_p;
		}
	}

	free(nodes);
	free(weights);
	free(f_wav_lmp);
	free(f_scal_lmp);
	free(flmp);
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
 * \param[in]  wav_lmp Wavelet kernel in FLAG space.
 * \param[in]  scal_lmp Scaling function kernel in FLAG space.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_p Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  P Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_p First wavelet scale to be used in radial space.
 * \retval none
 */
void flaglet_axisym_wav_synthesis_multires_real(double *f, const double *f_wav, const double *f_scal, double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p)
{
	complex double *f_wav_lmp, *f_scal_lmp;
	flaglet_axisym_allocate_f_wav_multires_lmp(&f_wav_lmp, &f_scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);

	int bandlimit_p, bandlimit_l, offset_lmp, offset, jl, jp;
	int J_l = s2let_j_max(L, B_l);
	int J_p = s2let_j_max(P, B_p);
	flag_core_analysis_real(f_scal_lmp, f_scal, R, L, P);

	offset_lmp = 0;
	offset = 0;
	for (jp = J_min_p; jp <= J_p; jp++){
		bandlimit_p = MIN(s2let_bandlimit(B_p, jp), P);
		for (jl = J_min_l; jl <= J_l; jl++){
			bandlimit_l = MIN(s2let_bandlimit(B_l, jl), L);
			//offset_lmp = jp * (J_l + 1) * L * L * P   +  jl * L * L * P;
			//offset = jp * (J_l + 1) * L * (2*L-1) * (P)   +  jl * L * (2*L-1) * (P);
			flag_core_analysis_real(f_wav_lmp + offset_lmp, f_wav + offset, R, bandlimit_l, bandlimit_p);
			offset_lmp += bandlimit_l * bandlimit_l * bandlimit_p ;
			offset += bandlimit_l * (2 * bandlimit_l - 1) * bandlimit_p;
		}
	}

	double *wav_lmp, *scal_lmp ;
	flaglet_axisym_allocate_wav_lmp(&wav_lmp, &scal_lmp, B_l, B_p, L, P);
	flaglet_axisym_wav_lmp(wav_lmp, scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);

	complex double *flmp;
	flag_core_allocate_flmn(&flmp, L, P);
	flaglet_axisym_wav_synthesis_multires_lmp(flmp, f_wav_lmp, f_scal_lmp, wav_lmp, scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);
	
	double *nodes = (double*)calloc(P, sizeof(double));
	double *weights = (double*)calloc(P, sizeof(double));
	assert(nodes != NULL);
	assert(weights != NULL);
	flag_spherlaguerre_sampling(nodes, weights, R, P);

	flag_core_synthesis_real(f, flmp, nodes, P, L, P);

	free(nodes);
	free(weights);
	free(wav_lmp);
	free(scal_lmp);
	free(f_wav_lmp);
	free(f_scal_lmp);
	free(flmp);
}

/*!
 * Perform wavelet transform in real space (from scratch, gives pixel space components).
 * Sampling scheme : MW sampling.
 * 3D spherical wavelets : analysis in real space.
 *
 * \param[out]  f_wav Wavelet transform (wavelet contribution in real space).
 * \param[out]  f_scal Wavelet transform (scaling contribution in real space).
 * \param[in]  f Input function (MW sampling).
 * \param[in]  wav_lmp Wavelet kernel in FLAG space.
 * \param[in]  scal_lmp Scaling function kernel in FLAG space.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_p Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  P Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_p First wavelet scale to be used in radial space.
 * \retval none
 */
void flaglet_axisym_wav_analysis(complex double *f_wav, complex double *f_scal, const complex double *f, double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p)
{
	complex double *flmp;
	flag_core_allocate_flmn(&flmp, L, P);
	flag_core_analysis(flmp, f, R, L, P);

	double *wav_lmp, *scal_lmp;
	complex double *f_wav_lmp, *f_scal_lmp;

	flaglet_axisym_allocate_wav_lmp(&wav_lmp, &scal_lmp, B_l, B_p, L, P);
	flaglet_axisym_wav_lmp(wav_lmp, scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);

	flaglet_axisym_allocate_f_wav_lmp(&f_wav_lmp, &f_scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);
	flaglet_axisym_wav_analysis_lmp(f_wav_lmp, f_scal_lmp, flmp, wav_lmp, scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);

	free(wav_lmp);
	free(scal_lmp);

	double *nodes = (double*)calloc(P, sizeof(double));
	double *weights = (double*)calloc(P, sizeof(double));
	assert(nodes != NULL);
	assert(weights != NULL);
	flag_spherlaguerre_sampling(nodes, weights, R, P);

	int offset_lmp, offset, jl, jp;
	int J_l = s2let_j_max(L, B_l);
	int J_p = s2let_j_max(P, B_p);
	flag_core_synthesis(f_scal, f_scal_lmp, nodes, P, L, P);

	offset_lmp = 0;
	offset = 0;
	for (jp = J_min_p; jp <= J_p; jp++){
		for (jl = J_min_l; jl <= J_l; jl++){
			//offset_lmp = jp * (J_l + 1) * L * L * P   +  jl * L * L * P;
			//offset = jp * (J_l + 1) * L * (2*L-1) * (P)  +  jl * L * (2*L-1) * (P);
			flag_core_synthesis(f_wav + offset, f_wav_lmp + offset_lmp, nodes, P, L, P);
			offset_lmp += L * L * P ;
			offset += L * (2 * L - 1) * P;
		}
	}

	free(nodes);
	free(weights);
	free(f_wav_lmp);
	free(f_scal_lmp);
	free(flmp);
}

/*!
 * Perform wavelet transform in real space (from scratch, gives pixel space components).
 * Sampling scheme : MW sampling.
 * 3D spherical wavelets : synthesis in real space.
 *
 * \param[out]  f Input function (MW sampling).
 * \param[in]  f_wav Wavelet transform (wavelet contribution in real space).
 * \param[in]  f_scal Wavelet transform (scaling contribution in real space).
 * \param[in]  wav_lmp Wavelet kernel in FLAG space.
 * \param[in]  scal_lmp Scaling function kernel in FLAG space.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_p Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  P Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_p First wavelet scale to be used in radial space.
 * \retval none
 */
void flaglet_axisym_wav_synthesis(complex double *f, const complex double *f_wav, const complex double *f_scal, double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p)
{
	complex double *f_wav_lmp, *f_scal_lmp;
	flaglet_axisym_allocate_f_wav_lmp(&f_wav_lmp, &f_scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);

	int offset_lmp, offset, jl, jp;
	int J_l = s2let_j_max(L, B_l);
	int J_p = s2let_j_max(P, B_p);
	flag_core_analysis(f_scal_lmp, f_scal, R, L, P);

	offset_lmp = 0;
	offset = 0;
	for (jp = J_min_p; jp <= J_p; jp++){
		for (jl = J_min_l; jl <= J_l; jl++){
			//offset_lmp = jp * (J_l + 1) * L * L * P   +  jl * L * L * P;
			//offset = jp * (J_l + 1) * L * (2*L-1) * (P)   +  jl * L * (2*L-1) * (P);
			flag_core_analysis(f_wav_lmp + offset_lmp, f_wav + offset, R, L, P);
			offset_lmp += L * L * P ;
			offset += L * (2 * L - 1) * P;
		}
	}

	double *wav_lmp, *scal_lmp ;
	flaglet_axisym_allocate_wav_lmp(&wav_lmp, &scal_lmp, B_l, B_p, L, P);
	flaglet_axisym_wav_lmp(wav_lmp, scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);

	complex double *flmp;
	flag_core_allocate_flmn(&flmp, L, P);
	flaglet_axisym_wav_synthesis_lmp(flmp, f_wav_lmp, f_scal_lmp, wav_lmp, scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);
	
	double *nodes = (double*)calloc(P, sizeof(double));
	double *weights = (double*)calloc(P, sizeof(double));
	assert(nodes != NULL);
	assert(weights != NULL);
	flag_spherlaguerre_sampling(nodes, weights, R, P);

	flag_core_synthesis(f, flmp, nodes, P, L, P);

	free(nodes);
	free(weights);
	free(wav_lmp);
	free(scal_lmp);
	free(f_wav_lmp);
	free(f_scal_lmp);
	free(flmp);
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
 * \param[in]  wav_lmp Wavelet kernel in FLAG space.
 * \param[in]  scal_lmp Scaling function kernel in FLAG space.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_p Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  P Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_p First wavelet scale to be used in radial space.
 * \retval none
 */
void flaglet_axisym_wav_analysis_real(double *f_wav, double *f_scal, const double *f, double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p)
{
	complex double *flmp;
	flag_core_allocate_flmn(&flmp, L, P);
	flag_core_analysis_real(flmp, f, R, L, P);

	double *wav_lmp, *scal_lmp;
	complex double *f_wav_lmp, *f_scal_lmp;

	flaglet_axisym_allocate_wav_lmp(&wav_lmp, &scal_lmp, B_l, B_p, L, P);
	flaglet_axisym_wav_lmp(wav_lmp, scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);

	flaglet_axisym_allocate_f_wav_lmp(&f_wav_lmp, &f_scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);
	flaglet_axisym_wav_analysis_lmp(f_wav_lmp, f_scal_lmp, flmp, wav_lmp, scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);

	free(wav_lmp);
	free(scal_lmp);

	double *nodes = (double*)calloc(P, sizeof(double));
	double *weights = (double*)calloc(P, sizeof(double));
	assert(nodes != NULL);
	assert(weights != NULL);
	flag_spherlaguerre_sampling(nodes, weights, R, P);

	int offset_lmp, offset, jl, jp;
	int J_l = s2let_j_max(L, B_l);
	int J_p = s2let_j_max(P, B_p);
	flag_core_synthesis_real(f_scal, f_scal_lmp, nodes, P, L, P);

	offset_lmp = 0;
	offset = 0;
	for (jp = J_min_p; jp <= J_p; jp++){
		for (jl = J_min_l; jl <= J_l; jl++){
			//offset_lmp = jp * (J_l + 1) * L * L * P   +  jl * L * L * P;
			//offset = jp * (J_l + 1) * L * (2*L-1) * (P)  +  jl * L * (2*L-1) * (P);
			flag_core_synthesis_real(f_wav + offset, f_wav_lmp + offset_lmp, nodes, P, L, P);
			offset_lmp += L * L * P ;
			offset += L * (2 * L - 1) * P;
		}
	}

	free(nodes);
	free(weights);
	free(f_wav_lmp);
	free(f_scal_lmp);
	free(flmp);
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
 * \param[in]  wav_lmp Wavelet kernel in FLAG space.
 * \param[in]  scal_lmp Scaling function kernel in FLAG space.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_p Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  P Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_p First wavelet scale to be used in radial space.
 * \retval none
 */
void flaglet_axisym_wav_synthesis_real(double *f, const double *f_wav, const double *f_scal, double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p)
{
	complex double *f_wav_lmp, *f_scal_lmp;
	flaglet_axisym_allocate_f_wav_lmp(&f_wav_lmp, &f_scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);

	int offset_lmp, offset, jl, jp;
	int J_l = s2let_j_max(L, B_l);
	int J_p = s2let_j_max(P, B_p);
	flag_core_analysis_real(f_scal_lmp, f_scal, R, L, P);

	offset_lmp = 0;
	offset = 0;
	for (jp = J_min_p; jp <= J_p; jp++){
		for (jl = J_min_l; jl <= J_l; jl++){
			//offset_lmp = jp * (J_l + 1) * L * L * P   +  jl * L * L * P;
			//offset = jp * (J_l + 1) * L * (2*L-1) * (P)   +  jl * L * (2*L-1) * (P);
			flag_core_analysis_real(f_wav_lmp + offset_lmp, f_wav + offset, R, L, P);
			offset_lmp += L * L * P ;
			offset += L * (2 * L - 1) * P;
		}
	}

	double *wav_lmp, *scal_lmp ;
	flaglet_axisym_allocate_wav_lmp(&wav_lmp, &scal_lmp, B_l, B_p, L, P);
	flaglet_axisym_wav_lmp(wav_lmp, scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);

	complex double *flmp;
	flag_core_allocate_flmn(&flmp, L, P);
	flaglet_axisym_wav_synthesis_lmp(flmp, f_wav_lmp, f_scal_lmp, wav_lmp, scal_lmp, B_l, B_p, L, P, J_min_l, J_min_p);
	
	double *nodes = (double*)calloc(P, sizeof(double));
	double *weights = (double*)calloc(P, sizeof(double));
	assert(nodes != NULL);
	assert(weights != NULL);
	flag_spherlaguerre_sampling(nodes, weights, R, P);

	flag_core_synthesis_real(f, flmp, nodes, P, L, P);

	free(nodes);
	free(weights);
	free(wav_lmp);
	free(scal_lmp);
	free(f_wav_lmp);
	free(f_scal_lmp);
	free(flmp);
}

/*!
 * Indice corresponding to a quintuplet (jl, jp, l, m, n) in the wavelets kernels.
 *
 * \param[in]  jl Angular scale indice.
 * \param[in]  jp Radial scale indice.
 * \param[in]  l Multipole indice.
 * \param[in]  m Order indice.
 * \param[in]  n Laguerre order indice.
 * \param[in]  J_l Maximum scale for angular harmonic space.
 * \param[in]  J_p Maximum scale for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  P Radial harmonic band-limit.
 * \retval Indice
 */
int jjlmp2ind(int jl, int jp, int l, int m, int n, int J_l, int J_p, int L, int P)
{
	return jp * (J_l + 1) * L * L * P   +  jl * L * L * P  +  n * L * L  +  l * l  + l +  m;
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
int lmp2ind(int l, int m, int n, int L)
{
	return n * L * L  +  l * l  +  l  +  m ;
}


