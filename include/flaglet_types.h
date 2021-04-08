#ifndef FLAGLET_TYPES
#define FLAGLET_TYPES

#include <s2let/s2let.h>
#include <so3/so3.h>

#ifdef __cplusplus
#include <complex>
#define FLAGLET_COMPLEX(TYPE) std::complex<TYPE>
extern "C" {
#else
#define FLAGLET_COMPLEX(TYPE) TYPE complex
#endif

typedef struct {

    int B_l;
    int L;
    int J_min_l;
    int N;

    int B_p;
    int P;
    int J_min_p;

    int spin;

    int upsample;

    int reality;

    double tau;

} flaglet_parameters_t;



static inline void fill_s2let_angular_parameters(s2let_parameters_t *s2let_parameters, const flaglet_parameters_t *parameters)
{
    s2let_parameters->upsample = parameters->upsample;
    s2let_parameters->L = parameters->L;
    s2let_parameters->N = parameters->N;
    s2let_parameters->spin = parameters->spin;
    s2let_parameters->B = parameters->B_l;
    s2let_parameters->J_min = parameters->J_min_l;
    s2let_parameters->sampling_scheme = S2LET_SAMPLING_MW;
    s2let_parameters->dl_method = SSHT_DL_RISBO;
    s2let_parameters->reality = parameters->reality;
}

static inline void fill_so3_angular_parameters(so3_parameters_t *so3_parameters, const flaglet_parameters_t *parameters)
{
    s2let_parameters_t s2let_parameters = {};
    fill_s2let_angular_parameters(&s2let_parameters, parameters);
    fill_so3_parameters(so3_parameters, &s2let_parameters);
}


#ifdef __cplusplus
}
#endif
#endif