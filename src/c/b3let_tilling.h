// B3LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef B3LET_TILLING
#define B3LET_TILLING

void b3let_axisym_allocate_tilling(double **kappa_ln, double **kappa0_ln, int B_l, int B_n, int L, int N);

void b3let_axisym_tilling(double *kappa_ln, double *kappa0_ln, int B_l, int B_n, int L, int N, int J_min_l, int J_min_n);

double b3let_axisym_check_identity(double *kappa_ln, double *kappa0_ln, int B_l, int B_n, int L, int N, int J_min_l, int J_min_n);

#endif