// FLAGLET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include "flaglet.h"
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <flag.h>
#include <s2let.h>
#include <so3.h>
#include <assert.h>


int flaglet_radial_bandlimit(int jp, const flaglet_parameters_t *parameters){
	return ceil(pow(parameters->B_p, jp+1));
}

int flaglet_angular_bandlimit(int jl, const flaglet_parameters_t *parameters){
	s2let_parameters_t s2let_parameters = {};
	fill_s2let_angular_parameters(&s2let_parameters, parameters);
	return s2let_bandlimit(jl, &s2let_parameters);
}



/*!
 * Allocates arrays for the kernels of the wavelets and the scaling functions (in FLAG space).
 */
void flaglet_allocate_wav_lmp(complex double **wav_lmp, double **scal_lmp, const flaglet_parameters_t *parameters)
{
	int J_l = flaglet_j_max(parameters->L, parameters->B_l);
	int J_p = flaglet_j_max(parameters->P, parameters->B_p);
	int L = parameters->L;
	int P = parameters->P;
	*wav_lmp = (complex double*)calloc( (J_l+1) * (J_p+1) * P * L * L, sizeof(complex double));
	*scal_lmp = (double*)calloc( L * P, sizeof(double));
}



/*!
 * Computes the kernels of the wavelets and the scaling functions (in FLAG space).
 */
void flaglet_wav_lmp(complex double *wav_lmp, double *scal_lmp, const flaglet_parameters_t *parameters)
{
	int L = parameters->L;
	int P = parameters->P;

	s2let_parameters_t s2let_parameters = {};
	fill_s2let_angular_parameters(&s2let_parameters, parameters);

	int jl, jp, l, p, m, indjjlmp, indlm;
	int J_l = flaglet_j_max(L, parameters->B_l);
	int J_p = flaglet_j_max(P, parameters->B_p);

	double wav0, scal0;

	double *kappa_ln, *kappa0_ln;
	flaglet_tiling_axisym_allocate(&kappa_ln, &kappa0_ln, parameters);

    complex double *s_elm;
    s2let_tiling_direction_allocate(&s_elm, &s2let_parameters);
    s2let_tiling_direction(s_elm, &s2let_parameters);

	flaglet_tiling_axisym(kappa_ln, kappa0_ln, parameters);

	for (jp = parameters->J_min_p; jp <= J_p; jp++){
		for (jl = parameters->J_min_l; jl <= J_l; jl++){
			for (p = 0; p < P; p++){
				for (l = 0; l < L; l++){
					wav0 = sqrt((2*l+1)/(8.0*PI*PI)) * kappa_ln[ jp*(J_l+1)*L*P  + jl*L*P + p*L + l ];
					for (m = -l; m <= l ; ++m){
						ssht_sampling_elm2ind(&indlm, l, m);
						indjjlmp = jp * (J_l + 1) * L * L * P   +  jl * L * L * P + p * L * L + indlm;
						wav_lmp[indjjlmp] =  wav0 * s_elm[indlm];
					}
				}
			}
		}
	}
	for (p = 0; p < P; p++){
		for (l = 0; l < L; l++){
			scal0 = sqrt( (2 * l + 1) / (4.0 * PI) ) * kappa0_ln[ p*L + l ];
			scal_lmp[p * L + l] = scal0;
		}
	}

	free(kappa_ln);
	free(kappa0_ln);
    free(s_elm);
}

double flaglet_tiling_wavelet_check_identity(complex double *wav_lmp, double *scal_lmp, const flaglet_parameters_t *parameters)
{
    int L = parameters->L;
    int P = parameters->P;
    int spin = parameters->spin;

    int indjjlmp, indlm, jl, jp, l, p, m;
	int J_l = flaglet_j_max(L, parameters->B_l);
	int J_p = flaglet_j_max(P, parameters->B_p);
    double error = 0.0; // maximum error for all el

    double *ident;
    ident = calloc(L*P, sizeof *ident);

	for (p = 0; p < P; p++){
		for (l = 0; l < L; l++){
			ident[p*L + l] += 4.0*PI/(2*l+1) * scal_lmp[p*L + l] * conj(scal_lmp[p*L + l]);
		}
	}

	for (jp = parameters->J_min_p; jp <= J_p; jp++){
		for (jl = parameters->J_min_l; jl <= J_l; jl++){
			for (p = 0; p < P; p++){
				for (l = 0; l < L; l++){
					for (m = -l; m <= l ; ++m){
						ssht_sampling_elm2ind(&indlm, l, m);
						indjjlmp = jp * (J_l + 1) * L * L * P   +  jl * L * L * P + p * L * L + indlm;
						ident[p*L+l] += 8.0*PI*PI/(2*l+1) * wav_lmp[indjjlmp] * conj(wav_lmp[indjjlmp]);
					}
				}
			}
		}
	}

	for (p = 0; p < P; p++){
		for (l = ABS(spin); l < L; l++){
			for (m = -l; m <= l ; m++){
        		error = MAX(error, fabs(ident[p*L+l] - 1.0));
		    	if( fabs(creal(ident[p*L+l])-1.0)>0.01)
					printf("(l,m,p) = (%i,%i,%i) - sq = (%f,%f)\n",l,m,p,creal(ident[p*L+l]),cimag(ident[p*L+l]));
		    }
		}
	}

    return error;
}

/*!
 * Perform multiresolution wavelet transform in FLAG space (from precomputed kernels, gives FLAG coefficients).
 * 3D spherical wavelets : analysis in FLAG-harmonic space.
 */
void flaglet_analysis_lmnp(complex double *f_wav_lmnp, complex double *f_scal_lmp, const complex double *flmp, const complex double *wav_lmp, const double *scal_lmp, const flaglet_parameters_t *parameters)
{
	int offset, jl, jp, l, m, n, p, lmn_size, lm_ind, ln_ind, indlmn, indjjlnp, indlmp;
	complex double psi;
	int L = parameters->L;
	int N = parameters->N;
	int P = parameters->P;
	int J_l = flaglet_j_max(L, parameters->B_l);
	int J_p = flaglet_j_max(P, parameters->B_p);
	so3_parameters_t so3_parameters = {};
	fill_so3_angular_parameters(&so3_parameters, parameters);
	int bandlimit_p = P;
	int bandlimit_l = L;
	int Nj = N;

	offset = 0;
	for (jp = parameters->J_min_p; jp <= J_p; jp++){
        if (!parameters->upsample)
        	bandlimit_p = MIN(flaglet_radial_bandlimit(jp, parameters), P);
		for (jl = parameters->J_min_l; jl <= J_l; jl++){
        	if (!parameters->upsample)
        	{
            	bandlimit_l = MIN(flaglet_angular_bandlimit(jl, parameters), L);
            	so3_parameters.L = bandlimit_l;
            	Nj = MIN(N,bandlimit_l);
            	Nj += (Nj+N)%2;
            	so3_parameters.N = Nj;
        	}
			lmn_size = so3_sampling_flmn_size(&so3_parameters);
			for (p = 0; p < bandlimit_p; p++){
	       		for (n = -Nj+1; n < Nj; n+=2){
					for (l = MAX(ABS(parameters->spin), ABS(n)); l < bandlimit_l; l++){
		       			ssht_sampling_elm2ind(&ln_ind, l, n);
		       			indjjlnp = jp * (J_l + 1) * L * L * P   +  jl * L * L * P + p * L * L + ln_ind;
						psi = 8*PI*PI/(2*l+1) * conj(wav_lmp[indjjlnp]);
						for (m = -l; m <= l ; m++){
	       					ssht_sampling_elm2ind(&lm_ind, l, m);
							indlmp = p * L * L + lm_ind;
							so3_sampling_elmn2ind(&indlmn, l, m, n, &so3_parameters);
							f_wav_lmnp[offset + p * lmn_size + indlmn] = flmp[indlmp] * psi ;
						}
					}
				}
			}
			offset += lmn_size * bandlimit_p;
		}
	}

	for (p = 0; p < P; p++){
		for (l = ABS(parameters->spin); l < L; l++){
			for (m = -l; m <= l ; m++){
				ssht_sampling_elm2ind(&lm_ind, l, m);
				indlmp = p * L * L  +  lm_ind;
				f_scal_lmp[indlmp] = flmp[indlmp] * sqrt((4.0*PI)/(2.0*l+1.0)) * scal_lmp[p*L+l] ;
			}
		}
	}
}



/*!
 * Perform multiresolution wavelet transform in FLAG space (from precomputed kernels, gives FLAG coefficients).
 * 3D spherical wavelets : synthesis in FLAG-harmonic space.
 */
void flaglet_synthesis_lmnp(complex double *flmp, const complex double *f_wav_lmnp, const complex double *f_scal_lmp, const complex double *wav_lmp, const double *scal_lmp, const flaglet_parameters_t *parameters)
{
	int offset, jl, jp, l, m, n, p, lmn_size, lm_ind, ln_ind, indlmn, indjjlnp, indlmp;
	complex double psi;
	int L = parameters->L;
	int N = parameters->N;
	int P = parameters->P;
	int J_l = flaglet_j_max(L, parameters->B_l);
	int J_p = flaglet_j_max(P, parameters->B_p);
	so3_parameters_t so3_parameters = {};
	fill_so3_angular_parameters(&so3_parameters, parameters);
	int bandlimit_p = P;
	int bandlimit_l = L;
	int Nj = N;

	offset = 0;
	for (jp = parameters->J_min_p; jp <= J_p; jp++){
        if (!parameters->upsample){
        	bandlimit_p = MIN(flaglet_radial_bandlimit(jp, parameters), P);
        }
		for (jl = parameters->J_min_l; jl <= J_l; jl++){
        	if (!parameters->upsample) {
            	bandlimit_l = MIN(flaglet_angular_bandlimit(jl, parameters), L);
            	so3_parameters.L = bandlimit_l;
            	Nj = MIN(N,bandlimit_l);
            	Nj += (Nj+N)%2;
            	so3_parameters.N = Nj;
        	}
	        so3_parameters.L0 = ceil(pow(parameters->B_l, jl-1));
			lmn_size = so3_sampling_flmn_size(&so3_parameters);
			for (p = 0; p < bandlimit_p; p++){
	       		for (n = -Nj+1; n < Nj; n+=2){
					for (l = MAX(ABS(parameters->spin), ABS(n)); l < bandlimit_l; l++){
		       			ssht_sampling_elm2ind(&ln_ind, l, n);
		       			indjjlnp = jp * (J_l + 1) * L * L * P   +  jl * L * L * P + p * L * L + ln_ind;
						psi = wav_lmp[indjjlnp];
						for (m = -l; m <= l ; m++){
	       					ssht_sampling_elm2ind(&lm_ind, l, m);
							indlmp = p * L * L + lm_ind;
							so3_sampling_elmn2ind(&indlmn, l, m, n, &so3_parameters);
							flmp[indlmp] += f_wav_lmnp[offset + p * lmn_size + indlmn] * psi ;
						}
					}
				}
			}
			offset += lmn_size * bandlimit_p;
		}
	}

	if (!parameters->upsample) {
		bandlimit_p = MIN(flaglet_radial_bandlimit(parameters->J_min_p-1, parameters), P);
		bandlimit_l = MIN(flaglet_angular_bandlimit(parameters->J_min_l-1, parameters), L);
	}
	for (p = 0; p < P; p++){
		for (l = ABS(parameters->spin); l < L; l++){
			for (m = -l; m <= l ; m++){
				ssht_sampling_elm2ind(&lm_ind, l, m);
				indlmp = p * L * L  +  lm_ind;
				flmp[indlmp] += f_scal_lmp[indlmp] * sqrt((4.0*PI)/(2.0*l+1.0)) * scal_lmp[p*L+l] ;
			}
		}
	}

}



/*!
 * Allocates 3D Wavelet transform in FLAG space.
 */
void flaglet_allocate_f_wav_lmnp(complex double **f_wav_lmnp, complex double **f_scal_lmp, const flaglet_parameters_t *parameters)
{
	int L = parameters->L;
	int P = parameters->P;
	int N = parameters->N;
	int J_l = flaglet_j_max(L, parameters->B_l);
	int J_p = flaglet_j_max(P, parameters->B_p);
	int jp, jl, total = 0;
	int Nj = N;
	so3_parameters_t so3_parameters = {};
	fill_so3_angular_parameters(&so3_parameters, parameters);
	int bandlimit_p = P;
	int bandlimit_l = L;
	for (jp = parameters->J_min_p; jp <= J_p; jp++){
		if (!parameters->upsample) {
            bandlimit_p = MIN(flaglet_radial_bandlimit(jp, parameters), P);
        }
		for (jl = parameters->J_min_l; jl <= J_l; jl++){
			if (!parameters->upsample) {
            	bandlimit_l = MIN(flaglet_angular_bandlimit(jl, parameters), L);
            	so3_parameters.L = bandlimit_l;
            	Nj = MIN(N,bandlimit_l);
            	Nj += (Nj+N)%2;
            	so3_parameters.N = Nj;
        	}
	        so3_parameters.L0 = ceil(pow(parameters->B_l, jl-1));
        	total += so3_sampling_flmn_size(&so3_parameters) * bandlimit_p ;
		}
	}
	*f_wav_lmnp = (complex double*)calloc( total, sizeof(complex double));
	*f_scal_lmp = (complex double*)calloc( L * L * P, sizeof(complex double));
}

/*!
 * Allocates multiresolution 3D Wavelet transform in real space.
 */
void flaglet_allocate_f_wav(complex double **f_wav, complex double **f_scal, const flaglet_parameters_t *parameters)
{
	int L = parameters->L;
	int P = parameters->P;
	int N = parameters->N;
	int J_l = flaglet_j_max(L, parameters->B_l);
	int J_p = flaglet_j_max(P, parameters->B_p);
	int jp, jl, total = 0;
	so3_parameters_t so3_parameters = {};
	fill_so3_angular_parameters(&so3_parameters, parameters);
	int Nj = MIN(N,L);
	int bandlimit_p = P;
	int bandlimit_l = L;
	for (jp = parameters->J_min_p; jp <= J_p; jp++){
		if (!parameters->upsample) {
            bandlimit_p = MIN(flaglet_radial_bandlimit(jp, parameters), P);
        }
		for (jl = parameters->J_min_l; jl <= J_l; jl++){
			if (!parameters->upsample) {
            	bandlimit_l = MIN(flaglet_angular_bandlimit(jl, parameters), L);
            	so3_parameters.L = bandlimit_l;
            	Nj = MIN(N,bandlimit_l);
            	Nj += (Nj+N)%2;
            	so3_parameters.N = Nj;
        	}
	        so3_parameters.L0 = ceil(pow(parameters->B_l, jl-1));
        	total += so3_sampling_f_size(&so3_parameters) * bandlimit_p ;
		}
	}
	*f_wav = (complex double*)calloc( total, sizeof(complex double));
	*f_scal = (complex double*)calloc( L * (2*L-1) * P, sizeof(complex double));
}



/*!
 * Perform wavelet transform in real space (from scratch, gives pixel space components).
 * Sampling scheme : MW sampling.
 * 3D spherical wavelets : analysis in real space.
 */
void flaglet_analysis(complex double *f_wav, complex double *f_scal, const complex double *f, const flaglet_parameters_t *parameters)
{
	double tau = parameters->tau;
	int L = parameters->L;
	int P = parameters->P;
	int N = parameters->N;
	int spin = parameters->spin;
	int J_l = flaglet_j_max(L, parameters->B_l);
	int J_p = flaglet_j_max(P, parameters->B_p);
	so3_parameters_t so3_parameters = {};
	fill_so3_angular_parameters(&so3_parameters, parameters);

	complex double *flmp, *f_wav_p;
	flag_core_allocate_flmn(&flmp, L, P);
	flag_core_analysis(flmp, f, L, tau, P, spin);

	complex double *wav_lmp;
	double *scal_lmp;
	flaglet_allocate_wav_lmp(&wav_lmp, &scal_lmp, parameters);
	flaglet_wav_lmp(wav_lmp, scal_lmp, parameters);

	complex double *f_wav_lmnp, *f_scal_lmp;
	flaglet_allocate_f_wav_lmnp(&f_wav_lmnp, &f_scal_lmp, parameters);
	flaglet_analysis_lmnp(f_wav_lmnp, f_scal_lmp, flmp, wav_lmp, scal_lmp, parameters);

	int offset_lmnp, offset, jl, jp, lmn_size, f_size, p;
	int Nj = MIN(N,L);

	double *nodes = (double*)calloc(P, sizeof(double));
	double *weights = (double*)calloc(P, sizeof(double));
	flag_spherlaguerre_sampling(nodes, weights, tau, P);

	flag_core_synthesis(f_scal, f_scal_lmp, nodes, P, L, tau, P, spin);

	int bandlimit_p = P;
	int bandlimit_l = L;

	offset_lmnp = 0;
	offset = 0;
	for (jp = parameters->J_min_p; jp <= J_p; jp++){
		if (!parameters->upsample) {
			bandlimit_p = MIN(flaglet_radial_bandlimit(jp, parameters), P);
			flag_spherlaguerre_sampling(nodes, weights, tau, bandlimit_p);
		}
		for (jl = parameters->J_min_l; jl <= J_l; jl++){
			if (!parameters->upsample) {
	            bandlimit_l = MIN(flaglet_angular_bandlimit(jl, parameters), L);
	            so3_parameters.L = bandlimit_l;
	            Nj = MIN(N,bandlimit_l);
	            Nj += (Nj+N)%2; // ensure N and Nj are both even or both odd
	            so3_parameters.N = Nj;
	        }
	        so3_parameters.L0 = ceil(pow(parameters->B_l, jl-1));
			lmn_size = so3_sampling_flmn_size(&so3_parameters);
			f_size = so3_sampling_f_size(&so3_parameters);
			f_wav_p = (complex double*)calloc( f_size * bandlimit_p, sizeof(complex double));
	        for (p = 0; p < bandlimit_p; p++){
		        so3_core_inverse_via_ssht(
		            f_wav_p + p * f_size,
		            f_wav_lmnp + offset_lmnp,
		            &so3_parameters
		        );
		        offset_lmnp += lmn_size;
		    }
	        flag_spherlaguerre_mapped_synthesis(f_wav + offset, f_wav_p, nodes, bandlimit_p, tau, bandlimit_p, f_size);
			offset += f_size * bandlimit_p;
	        free(f_wav_p);
		}
	}

	free(nodes);
	free(weights);
	free(f_wav_lmnp);
	free(f_scal_lmp);
	free(flmp);
}

/*!
 * Perform wavelet transform in real space (from scratch, gives pixel space components).
 * Sampling scheme : MW sampling.
 * 3D spherical wavelets : synthesis in real space.
 */
void flaglet_synthesis(complex double *f, const complex double *f_wav, const complex double *f_scal, const flaglet_parameters_t *parameters)
{
	complex double *f_wav_lmnp, *f_scal_lmp, *f_wav_p;
	flaglet_allocate_f_wav_lmnp(&f_wav_lmnp, &f_scal_lmp, parameters);

	int Nj, offset_lmnp, offset, jl, jp, lmn_size, f_size, p;
	double tau = parameters->tau;
	int L = parameters->L;
	int P = parameters->P;
	int N = parameters->N;
	int spin = parameters->spin;
	int J_l = flaglet_j_max(L, parameters->B_l);
	int J_p = flaglet_j_max(P, parameters->B_p);
	so3_parameters_t so3_parameters = {};
	fill_so3_angular_parameters(&so3_parameters, parameters);

	double *nodes = (double*)calloc(P, sizeof(double));
	double *weights = (double*)calloc(P, sizeof(double));

	int bandlimit_p = P;
	int bandlimit_l = L;
	flag_spherlaguerre_sampling(nodes, weights, tau, P);

	flag_core_analysis(f_scal_lmp, f_scal, L, tau, P, spin);

	offset_lmnp = 0;
	offset = 0;
	for (jp = parameters->J_min_p; jp <= J_p; jp++){
		if (!parameters->upsample) {
			bandlimit_p = MIN(flaglet_radial_bandlimit(jp, parameters), P);
			flag_spherlaguerre_sampling(nodes, weights, tau, bandlimit_p);
		}
		for (jl = parameters->J_min_l; jl <= J_l; jl++){
			 if (!parameters->upsample){
	            bandlimit_l = MIN(flaglet_angular_bandlimit(jl, parameters), L);
	            so3_parameters.L = bandlimit_l;
	            Nj = MIN(N,bandlimit_l);
	            Nj += (Nj+N)%2; // ensure N and Nj are both even or both odd
	            so3_parameters.N = Nj;
	        }
	        so3_parameters.L0 = ceil(pow(parameters->B_l, jl-1));
			lmn_size = so3_sampling_flmn_size(&so3_parameters);
			f_size = so3_sampling_f_size(&so3_parameters);

			f_wav_p = (complex double*)calloc( f_size * bandlimit_p, sizeof(complex double));
			flag_spherlaguerre_mapped_analysis(f_wav_p, f_wav + offset, nodes, weights, tau, bandlimit_p, f_size);
			for (p = 0; p < bandlimit_p; p++){
		        so3_core_forward_via_ssht(
		            f_wav_lmnp + offset_lmnp,
		            f_wav_p + p * f_size,
		            &so3_parameters
		        );
		        offset_lmnp += lmn_size;
		    }
		    offset += f_size * bandlimit_p;
	        free(f_wav_p);
		}
	}

	complex double *wav_lmp;
	double *scal_lmp;
	flaglet_allocate_wav_lmp(&wav_lmp, &scal_lmp, parameters);
	flaglet_wav_lmp(wav_lmp, scal_lmp, parameters);

	complex double *flmp;
	flag_core_allocate_flmn(&flmp, L, P);
	flaglet_synthesis_lmnp(flmp, f_wav_lmnp, f_scal_lmp, wav_lmp, scal_lmp, parameters);

	flag_spherlaguerre_sampling(nodes, weights, tau, P);
	flag_core_synthesis(f, flmp, nodes, P, L, tau, P, spin);

	free(nodes);
	free(weights);
	free(wav_lmp);
	free(scal_lmp);
	free(f_wav_lmnp);
	free(f_scal_lmp);
	free(flmp);
}


