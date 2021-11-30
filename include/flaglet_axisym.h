// FLAGLET package
// Copyright (C) 2021
// Boris Leistedt & Jason McEwen & Matthew Price

#ifndef FLAGLET_AXISYM
#define FLAGLET_AXISYM

#include <complex.h> 
#ifdef __cplusplus
extern "C" {
#endif

int s2let_bandlimit_wrapper(int j, int J_min, int B, int L);
int s2let_j_max_wrapper(int L, int B);

void flaglet_axisym_allocate_f_wav_lmp(FLAGLET_COMPLEX(double) **f_wav_lmp, FLAGLET_COMPLEX(double) **f_scal_lmp, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);

void flaglet_axisym_allocate_f_wav_multires_lmp(FLAGLET_COMPLEX(double) **f_wav_lmp, FLAGLET_COMPLEX(double) **f_scal_lmp, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);

void flaglet_axisym_allocate_f_wav(FLAGLET_COMPLEX(double) **f_wav, FLAGLET_COMPLEX(double) **f_scal, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);
void flaglet_axisym_allocate_f_wav_real(double **f_wav, double **f_scal, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);

void flaglet_axisym_allocate_f_wav_multires(FLAGLET_COMPLEX(double) **f_wav, FLAGLET_COMPLEX(double) **f_scal, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);
void flaglet_axisym_allocate_f_wav_multires_real(double **f_wav, double **f_scal, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);

void flaglet_axisym_allocate_wav_lmp(double **wav_lmp, double **scal_lmp, int B_l, int B_p, int L, int P);
void flaglet_axisym_wav_lmp(double *wav_lmp, double *scal_lmp, int B_l, int B_p, int J_min_l, int J_min_p, int L, int P);

void flaglet_axisym_wav_analysis_lmp(FLAGLET_COMPLEX(double) *f_wav_lmp, FLAGLET_COMPLEX(double) *f_scal_lmp, const FLAGLET_COMPLEX(double) *flmp, const double *wav_lmp, const double *scal_lmp, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);
void flaglet_axisym_wav_synthesis_lmp(FLAGLET_COMPLEX(double) *flmp, const FLAGLET_COMPLEX(double) *f_wav_lmp, const FLAGLET_COMPLEX(double) *f_scal_lmp, const double *wav_lmp, const double *scal_lmp, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);
 	
void flaglet_axisym_wav_analysis_multires_lmp(FLAGLET_COMPLEX(double) *f_wav_lmp, FLAGLET_COMPLEX(double) *f_scal_lmp, const FLAGLET_COMPLEX(double) *flmp, const double *wav_lmp, const double *scal_lmp, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);
void flaglet_axisym_wav_synthesis_multires_lmp(FLAGLET_COMPLEX(double) *flmp, const FLAGLET_COMPLEX(double) *f_wav_lmp, const FLAGLET_COMPLEX(double) *f_scal_lmp, const double *wav_lmp, const double *scal_lmp, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);

void flaglet_axisym_wav_analysis(FLAGLET_COMPLEX(double) *f_wav, FLAGLET_COMPLEX(double) *f_scal, const FLAGLET_COMPLEX(double) *f, double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);
void flaglet_axisym_wav_synthesis(FLAGLET_COMPLEX(double) *f, const FLAGLET_COMPLEX(double) *f_wav, const FLAGLET_COMPLEX(double) *f_scal, double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);
void flaglet_axisym_wav_analysis_real(double *f_wav, double *f_scal, const double *f, double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);
void flaglet_axisym_wav_synthesis_real(double *f, const double *f_wav, const double *f_scal, double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);

void flaglet_axisym_wav_analysis_multires(FLAGLET_COMPLEX(double) *f_wav, FLAGLET_COMPLEX(double) *f_scal, const FLAGLET_COMPLEX(double) *f, double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);
void flaglet_axisym_wav_synthesis_multires(FLAGLET_COMPLEX(double) *f, const FLAGLET_COMPLEX(double) *f_wav, const FLAGLET_COMPLEX(double) *f_scal, double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);
void flaglet_axisym_wav_analysis_multires_real(double *f_wav, double *f_scal, const double *f, double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);
void flaglet_axisym_wav_synthesis_multires_real(double *f, const double *f_wav, const double *f_scal, double R, int B_l, int B_p, int L, int P, int J_min_l, int J_min_p);

int jjlmp2ind(int jl, int jn, int l, int m, int n, int J_l, int J_p, int L, int P);
int lmp2ind(int l, int m, int n, int L);

#ifdef __cplusplus
}
#endif
#endif
