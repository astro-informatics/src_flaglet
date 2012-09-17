// FLAGLET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef FLAGLET_AXISYM
#define FLAGLET_AXISYM

#include <complex.h> 

void flaglet_axisym_allocate_f_wav_lmp(complex double **f_wav_lmp, complex double **f_scal_lmp, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);

void flaglet_axisym_allocate_f_wav_multires_lmp(complex double **f_wav_lmp, complex double **f_scal_lmp, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);

void flaglet_axisym_allocate_f_wav(complex double **f_wav, complex double **f_scal, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);
void flaglet_axisym_allocate_f_wav_real(double **f_wav, double **f_scal, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);

void flaglet_axisym_allocate_f_wav_multires(complex double **f_wav, complex double **f_scal, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);
void flaglet_axisym_allocate_f_wav_multires_real(double **f_wav, double **f_scal, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);

void flaglet_axisym_allocate_wav_lmp(double **wav_lmp, double **scal_lmp, int B_l, int B_p, int L, int P);
void flaglet_axisym_wav_lmp(double *wav_lmp, double *scal_lmp, int B_l, int B_p, int J_min_l, int J_min_p, int L, int P);

void flaglet_axisym_wav_analysis_lmp(complex double *f_wav_lmp, complex double *f_scal_lmp, const complex double *flmp, const double *wav_lmp, const double *scal_lmp, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);
void flaglet_axisym_wav_synthesis_lmp(complex double *flmp, const complex double *f_wav_lmp, const complex double *f_scal_lmp, const double *wav_lmp, const double *scal_lmp, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);
 	
void flaglet_axisym_wav_analysis_multires_lmp(complex double *f_wav_lmp, complex double *f_scal_lmp, const complex double *flmp, const double *wav_lmp, const double *scal_lmp, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);
void flaglet_axisym_wav_synthesis_multires_lmp(complex double *flmp, const complex double *f_wav_lmp, const complex double *f_scal_lmp, const double *wav_lmp, const double *scal_lmp, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);

void flaglet_axisym_wav_analysis(complex double *f_wav, complex double *f_scal, const complex double *f, double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);
void flaglet_axisym_wav_synthesis(complex double *f, const complex double *f_wav, const complex double *f_scal, double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);
void flaglet_axisym_wav_analysis_real(double *f_wav, double *f_scal, const double *f, double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);
void flaglet_axisym_wav_synthesis_real(double *f, const double *f_wav, const double *f_scal, double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);

void flaglet_axisym_wav_analysis_multires(complex double *f_wav, complex double *f_scal, const complex double *f, double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);
void flaglet_axisym_wav_synthesis_multires(complex double *f, const complex double *f_wav, const complex double *f_scal, double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);
void flaglet_axisym_wav_analysis_multires_real(double *f_wav, double *f_scal, const double *f, double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);
void flaglet_axisym_wav_synthesis_multires_real(double *f, const double *f_wav, const double *f_scal, double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);

int jjlmp2ind(int jl, int jn, int l, int m, int n, int J_l, int J_p, int L, int P);
int lmp2ind(int l, int m, int n, int L);

#endif
