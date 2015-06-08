// FLAGLET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef FLAGLET_TRANSFORM
#define FLAGLET_TRANSFORM

#include <complex.h> 

int flaglet_radial_bandlimit(int jp, const flaglet_parameters_t *parameters);

int flaglet_angular_bandlimit(int jl, const flaglet_parameters_t *parameters);

void flaglet_allocate_f_wav_lmnp(complex double **f_wav_lmp, complex double **f_scal_lmp, const flaglet_parameters_t *parameters);

void flaglet_allocate_wav_lmp(complex double **wav_lmp, double **scal_lmp, const flaglet_parameters_t *parameters);

void flaglet_wav_lmp(complex double *wav_lmp, double *scal_lmp, const flaglet_parameters_t *parameters);

void flaglet_allocate_f_wav(complex double **f_wav, complex double **f_scal, const flaglet_parameters_t *parameters);

void flaglet_analysis_lmnp(complex double *f_wav_lmnp, complex double *f_scal_lmp, const complex double *flmp, const complex double *wav_lmp, const double *scal_lmp, const flaglet_parameters_t *parameters);
void flaglet_synthesis_lmnp(complex double *flmp, const complex double *f_wav_lmnp, const complex double *f_scal_lmp, const complex double *wav_lmp, const double *scal_lmp, const flaglet_parameters_t *parameters);
 	
void flaglet_analysis(complex double *f_wav, complex double *f_scal, const complex double *f, const flaglet_parameters_t *parameters);
void flaglet_synthesis(complex double *f, const complex double *f_wav, const complex double *f_scal, const flaglet_parameters_t *parameters);



#endif
