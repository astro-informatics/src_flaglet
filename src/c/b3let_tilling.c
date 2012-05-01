// B3LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "b3let.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))

/*!
 * Allocates tilling in FLAG - harmonic space.
 *
 * \param[out]  kappa_lp Kernel functions for the wavelets.
 * \param[out]  kappa0_lp Kernel for the scaling function.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_p Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  P Radial harmonic band-limit.
 * \retval none
 */
void b3let_axisym_allocate_tilling(double **kappa_lp, double **kappa0_lp, int B_l, int B_p, int L, int P)
{
	int J_l = s2let_j_max(L, B_l);
	int J_p = s2let_j_max(P, B_p);
	*kappa_lp = (double*)calloc( (J_p+1) * (J_l+1) * L * P, sizeof(double));
	*kappa0_lp = (double*)calloc( L * P, sizeof(double));
}

/*!
 * Allocates tilling in FLAG - harmonic space.
 *
 * \param[out]  kappa_lp Kernel functions for the wavelets.
 * \param[out]  kappa0_lp Kernel for the scaling function.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_p Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  P Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_p First wavelet scale to be used in radial space.
 * \retval none
 */
void b3let_axisym_tilling(double *kappa_lp, double *kappa0_lp, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p)
{
	int l, n, jl, jp;
	int J_l = s2let_j_max(L, B_l);
	int J_p = s2let_j_max(P, B_p);
	double temp;

	double *phi2_l = (double*)calloc((J_l+2) * L, sizeof(double));
	double *phi2_n = (double*)calloc((J_p+2) * P, sizeof(double));

	s2let_tilling_phi2(phi2_l, B_l, L, J_min_l);
	s2let_tilling_phi2(phi2_n, B_p, P, J_min_p);

	int el_max = ceil(pow(B_l,J_min_l))+1;
	int en_max = ceil(pow(B_p,J_min_p))+1;

	for (n = 0; n < en_max; n++){
		for (l = 0; l < el_max; l++){
			temp = sqrt(
				phi2_l[l+J_min_l*L] * phi2_n[n+(J_min_p+1)*P]
				+ phi2_l[l+(J_min_l+1)*L] * phi2_n[n+J_min_p*P]
				- phi2_l[l+J_min_l*L] * phi2_n[n+J_min_p*P]
			);
			if( isnan(temp) || isinf(temp) )
				kappa0_lp[n * L + l] = 0.0;
			else
				kappa0_lp[n * L + l] = temp;
		}
	}
	for (n = en_max; n < P; n++){
		for (l = 0; l < L; l++){
			temp = sqrt(phi2_l[l+J_min_l*L]);
			if( isnan(temp) || isinf(temp) )
				kappa0_lp[n * L + l] = 0.0;
			else
				kappa0_lp[n * L + l] = temp;
		}
	}
	for (n = 0; n < P; n++){
		for (l = el_max; l < L; l++){
			temp = sqrt(phi2_n[n+J_min_p*P]);
			if( isnan(temp) || isinf(temp) )
				kappa0_lp[n * L + l] = 0.0;
			else
				kappa0_lp[n * L + l] = temp;
		}
	}

	/*for (n = 0; n < P; n++){
		for (l = 0; l < L; l++){
			printf("  %f  ",kappa0_lp[n * L + l]);
		}
		printf("\n");
	}*/

	for (jp = J_min_p; jp <= J_p; jp++){
		for (jl = J_min_l; jl <= J_l; jl++){
			//printf("- j_l = %i - j_p = %i -\n",jl,jp);
			for (n = 0; n < P; n++){
				for (l = 0; l < L; l++){
					temp =
						sqrt(phi2_l[l+(jl+1)*L] - phi2_l[l+jl*L])
						* sqrt(phi2_n[n+(jp+1)*P] - phi2_n[n+jp*P]);
					if( isnan(temp) || isinf(temp) ){
						kappa_lp[ jp*(J_l+1)*L*P  + jl*L*P + n*L + l ] = 0.0;
					}else
						kappa_lp[ jp*(J_l+1)*L*P  + jl*L*P + n*L + l ] = temp;
					
					//printf(" %f ",kappa_lp[ jp*(J_l+1)*L*P  + jl*L*P + n*L + l ]);
				}
				//printf("\n");
			}
		}
	}

	free(phi2_l);
	free(phi2_n);
}

/*!
 * Checks exactness of the FLAG - harmonic tilling.
 *
 * \param[in]  kappa_lp Kernel functions for the wavelets.
 * \param[in]  kappa0_lp Kernel for the scaling function.
 * \param[in]  B_l Wavelet parameter for angular harmonic space.
 * \param[in]  B_p Wavelet parameter for radial harmonic space.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  P Radial harmonic band-limit.
 * \param[in]  J_min_l First wavelet scale to be used in angular space.
 * \param[in]  J_min_p First wavelet scale to be used in radial space.
 * \retval Achieved accuracy (should be lower than e-12).
 */
double b3let_axisym_check_identity(double *kappa_lp, double *kappa0_lp, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p)
{
	int jl, jp, l, n;
	int J_l = s2let_j_max(L, B_l);
	int J_p = s2let_j_max(P, B_p);
		
	double *ident;
	ident = (double*)calloc(L * P, sizeof(double));

	//printf("Tilling identity: \n");
	for (n = 0; n < P; n++){
		for (l = 0; l < L; l++){
			ident[l+n*L] = pow(kappa0_lp[n * L + l], 2.0);
		}
	}

	double sum = 0;
	for (n = 0; n < P; n++){
		for (l = 0; l < L; l++){
			for (jp = J_min_p; jp <= J_p; jp++){
				for (jl = J_min_l; jl <= J_l; jl++){
					//printf(" %f ",kappa_lp[ jp*(J_l+1)*L*P  + jl*L*P + n*L + l ]);
					ident[l+n*L] += pow(kappa_lp[ jp*(J_l+1)*L*P  + jl*L*P + n*L + l ], 2.0);
				}
			}
			//printf(" %2.2f ", ident[l+n*L]);
			sum += ident[l+n*L] - 1.000;
		}
		//printf("\n");
	}

	free(ident);
	return sum;
}



