// B3LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef FLAGLET_TILING
#define FLAGLET_TILING

#include <complex.h> 
#ifdef __cplusplus
extern "C" {
#endif
	
int flaglet_j_max(int L, int B);

void flaglet_tiling_axisym_allocate(double **kappa_lp, double **kappa0_lp, const flaglet_parameters_t *parameters);

void flaglet_tiling_axisym(double *kappa_lp, double *kappa0_lp, const flaglet_parameters_t *parameters);

double flaglet_tiling_axisym_check_identity(double *kappa_lp, double *kappa0_lp, const flaglet_parameters_t *parameters);

#ifdef __cplusplus
}
#endif
#endif
