// B3LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef B3LET_AXISYM
#define B3LET_AXISYM

void b3let_allocate_wav_lmn(double **wav_lmn, double **scal_lmn, int B_l, int B_n, int L, int N);
void b3let_wav_lmn(double *wav_lmn, double *scal_lmn, int B_l, int B_n, int J_min_l, int J_min_n, int L, int N);

void b3let_allocate_f_wav_lmn(complex double **f_wav_lmn, complex double **f_scal_lmn, int B_l, int B_n, int L, int N);
void b3let_wav_analysis_lmn(complex double *f_wav_lmn, complex double *f_scal_lmn, const complex double *flmn, const double *wav_lmn, const double *scal_lmn, int B_l, int B_n, int L, int N, int J_min_l, int J_min_n);
void b3let_wav_synthesis_lmn(complex double *flmn, const complex double *f_wav_lmn, const complex double *f_scal_lmn, const double *wav_lmn, const double *scal_lmn, int B_l, int B_n, int L, int N, int J_min_l, int J_min_n);
 	
void b3let_wav_analysis(complex double *f_wav, complex double *f_scal, const complex double *f, int B_l, int B_n, int L, int N, int J_min_l, int J_min_n);
void b3let_wav_synthesis(complex double *f, const complex double *f_wav, const complex double *f_scal, int B_l, int B_n, int L, int N, int J_min_l, int J_min_n);

int jjlmn2ind(int jl, int jn, int l, int m, int n, int J_l, int J_n, int L, int N);
int lmn2ind(int l, int m, int n, int L);

#endif