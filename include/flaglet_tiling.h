// B3LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef FLAGLET_TILING
#define FLAGLET_TILING

#include <complex.h> 

void flaglet_tiling_axisym_allocate(double **kappa_lp, double **kappa0_lp, int B_l, int B_p, int L, int P);

void flaglet_tiling_axisym(double *kappa_lp, double *kappa0_lp, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);

double flaglet_tiling_axisym_check_identity(double *kappa_lp, double *kappa0_lp, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);

#endif
