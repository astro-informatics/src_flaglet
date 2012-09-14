// B3LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef FLAGLET_TILLING
#define FLAGLET_TILLING

void flaglet_axisym_allocate_tilling(double **kappa_lp, double **kappa0_lp, int B_l, int B_p, int L, int P);

void flaglet_axisym_tilling(double *kappa_lp, double *kappa0_lp, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);

double flaglet_axisym_check_identity(double *kappa_lp, double *kappa0_lp, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);

#endif