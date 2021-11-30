// FLAGLET package
// Copyright (C) 2021
// Boris Leistedt & Jason McEwen & Matthew Price

#ifndef FLAGLET_TRANSFORM
#define FLAGLET_TRANSFORM

#include <complex.h> 
#ifdef __cplusplus
extern "C" {
#endif

int flaglet_radial_bandlimit(int jp, const flaglet_parameters_t *parameters);
int flaglet_angular_bandlimit(int jl, const flaglet_parameters_t *parameters);

int flaglet_f_size(const flaglet_parameters_t *parameters);
int flaglet_scal_size(const flaglet_parameters_t *parameters);
int flaglet_wav_size(const flaglet_parameters_t *parameters);

void flaglet_allocate_f_wav_lmnp(FLAGLET_COMPLEX(double) **f_wav_lmp, FLAGLET_COMPLEX(double) **f_scal_lmp, const flaglet_parameters_t *parameters);
void flaglet_allocate_wav_lmp(FLAGLET_COMPLEX(double) **wav_lmp, double **scal_lmp, const flaglet_parameters_t *parameters);
void flaglet_wav_lmp(FLAGLET_COMPLEX(double) *wav_lmp, double *scal_lmp, const flaglet_parameters_t *parameters);
void flaglet_allocate_f_wav(FLAGLET_COMPLEX(double) **f_wav, FLAGLET_COMPLEX(double) **f_scal, const flaglet_parameters_t *parameters);
void flaglet_analysis_lmnp(FLAGLET_COMPLEX(double) *f_wav_lmnp, FLAGLET_COMPLEX(double) *f_scal_lmp, const FLAGLET_COMPLEX(double) *flmp, const FLAGLET_COMPLEX(double) *wav_lmp, const double *scal_lmp, const flaglet_parameters_t *parameters);
void flaglet_analysis_lmnp_adjoint(FLAGLET_COMPLEX(double) *flmp, const FLAGLET_COMPLEX(double) *f_wav_lmnp, const FLAGLET_COMPLEX(double) *f_scal_lmp, const FLAGLET_COMPLEX(double) *wav_lmp, const double *scal_lmp, const flaglet_parameters_t *parameters);
void flaglet_synthesis_lmnp(FLAGLET_COMPLEX(double) *flmp, const FLAGLET_COMPLEX(double) *f_wav_lmnp, const FLAGLET_COMPLEX(double) *f_scal_lmp, const FLAGLET_COMPLEX(double) *wav_lmp, const double *scal_lmp, const flaglet_parameters_t *parameters);
void flaglet_synthesis_lmnp_adjoint(FLAGLET_COMPLEX(double) *f_wav_lmnp, FLAGLET_COMPLEX(double) *f_scal_lmp, const FLAGLET_COMPLEX(double) *flmp, const FLAGLET_COMPLEX(double) *wav_lmp, const double *scal_lmp, const flaglet_parameters_t *parameters);
void flaglet_analysis(FLAGLET_COMPLEX(double) *f_wav, FLAGLET_COMPLEX(double) *f_scal, const FLAGLET_COMPLEX(double) *f, const flaglet_parameters_t *parameters);
void flaglet_analysis_adjoint(FLAGLET_COMPLEX(double) *f, const FLAGLET_COMPLEX(double) *f_wav, const FLAGLET_COMPLEX(double) *f_scal, const flaglet_parameters_t *parameters);
void flaglet_synthesis(FLAGLET_COMPLEX(double) *f, const FLAGLET_COMPLEX(double) *f_wav, const FLAGLET_COMPLEX(double) *f_scal, const flaglet_parameters_t *parameters);
void flaglet_synthesis_adjoint(FLAGLET_COMPLEX(double) *f_wav, FLAGLET_COMPLEX(double) *f_scal, const FLAGLET_COMPLEX(double) *f, const flaglet_parameters_t *parameters);


#ifdef __cplusplus
}
#endif
#endif
