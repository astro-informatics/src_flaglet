// B3LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "b3let.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))

void b3let_allocate_tilling(double **kappa_ln, double **kappa0_ln, int B_l, int B_n, int L, int N)
{
	int J_l = s2let_j_max(L, B_l);
	int J_n = s2let_j_max(N, B_n);
	*kappa_ln = (double*)calloc( (J_n+1) * (J_l+1) * L * N, sizeof(double));
	*kappa0_ln = (double*)calloc( L * N, sizeof(double));
}

void b3let_tilling(double *kappa_ln, double *kappa0_ln, int B_l, int B_n, int L, int N, int J_min_l, int J_min_n)
{
	int l, m, n, jl, jn;
	int J_l = s2let_j_max(L, B_l);
	int J_n = s2let_j_max(N, B_n);

	double *phi2_l = (double*)calloc((J_l+2) * L, sizeof(double));
	double *phi2_n = (double*)calloc((J_n+2) * N, sizeof(double));

	s2let_tilling_phi2(phi2_l, B_l, L, J_min_l);
	s2let_tilling_phi2(phi2_n, B_n, N, J_min_n);

	int el_max = ceil(pow(B_l,J_min_l))+1;
	int en_max = ceil(pow(B_n,J_min_n))+1;

	for (n = 0; n < en_max; n++){
		for (l = 0; l < el_max; l++){
			kappa0_ln[n * L + l] = sqrt(
				phi2_l[l+J_min_l*L] * phi2_n[n+(J_min_n+1)*N]
				+ phi2_l[l+(J_min_l+1)*L] * phi2_n[n+J_min_n*N]
				- phi2_l[l+J_min_l*L] * phi2_n[n+J_min_n*N]
			);
		}
	}
	for (n = en_max; n < N; n++){
		for (l = 0; l < L; l++){
			kappa0_ln[n * L + l] = sqrt(phi2_l[l+J_min_l*L]);
		}
	}
	for (n = 0; n < N; n++){
		for (l = el_max; l < L; l++){
			kappa0_ln[n * L + l] = sqrt(phi2_n[n+J_min_n*N]);
		}
	}

	/*for (n = 0; n < N; n++){
		for (l = 0; l < L; l++){
			printf("  %f  ",kappa0_ln[n * L + l]);
		}
		printf("\n");
	}*/

	for (jn = J_min_n; jn <= J_n; jn++){
		for (jl = J_min_l; jl <= J_l; jl++){
			//printf("- j_l = %i - j_n = %i -\n",jl,jn);
			for (n = 0; n < N; n++){
				for (l = 0; l < L; l++){
					kappa_ln[ jn*(J_l+1)*L*N  + jl*L*N + n*L + l ] =
						sqrt(phi2_l[l+(jl+1)*L] - phi2_l[l+jl*L])
						* sqrt(phi2_n[n+(jn+1)*N] - phi2_n[n+jn*N]);
					//printf(" %f ",kappa_ln[ jn*(J_l+1)*L*N  + jl*L*N + n*L + l ]);
				}
				//printf("\n");
			}
		}
	}

	free(phi2_l);
	free(phi2_n);
}

double b3let_check_identity(double *kappa_ln, double *kappa0_ln, int B_l, int B_n, int L, int N, int J_min_l, int J_min_n)
{
	int jl, jn, l, n;
	int J_l = s2let_j_max(L, B_l);
	int J_n = s2let_j_max(N, B_n);
		
	double *ident;
	ident = (double*)calloc(L * N, sizeof(double));

	//printf("Tilling identity: \n");
	for (n = 0; n < N; n++){
		for (l = 0; l < L; l++){
			ident[l+n*L] = pow(kappa0_ln[n * L + l], 2.0);
		}
	}

	double sum = 0;
	for (n = 0; n < N; n++){
		for (l = 0; l < L; l++){
			for (jn = J_min_n; jn <= J_n; jn++){
				for (jl = J_min_l; jl <= J_l; jl++){
					//printf(" %f ",kappa_ln[ jn*(J_l+1)*L*N  + jl*L*N + n*L + l ]);
					ident[l+n*L] += pow(kappa_ln[ jn*(J_l+1)*L*N  + jl*L*N + n*L + l ], 2.0);
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



