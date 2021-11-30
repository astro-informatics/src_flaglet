// FLAGLET package
// Copyright (C) 2021
// Boris Leistedt & Jason McEwen & Matthew Price

void flaglet_analysis_adjoint(complex double *f, const complex double *f_wav, const complex double *f_scal, const flaglet_parameters_t *parameters)
{
	double tau = parameters->tau;
	int L = parameters->L;
	int P = parameters->P;
	int N = parameters->N;
	int spin = parameters->spin;
	int J_l = flaglet_j_max(L, parameters->B_l);
	int J_p = flaglet_j_max(P, parameters->B_p);
	so3_parameters_t so3_parameters = {};
	fill_so3_angular_parameters(&so3_parameters, parameters);

	complex double *f_wav_lmnp, *f_scal_lmp;
	flaglet_allocate_f_wav_lmnp(&f_wav_lmnp, &f_scal_lmp, parameters);

	complex double *wav_lmp;
	double *scal_lmp;
	flaglet_allocate_wav_lmp(&wav_lmp, &scal_lmp, parameters);

	complex double *flmp, *f_wav_p;
	flag_core_allocate_flmn(&flmp, L, P);

	int offset_lmnp, offset, jl, jp, lmn_size, f_size, p;
	int Nj = MIN(N,L);

	double *nodes = (double*)calloc(P, sizeof(double));
	double *weights = (double*)calloc(P, sizeof(double));
	flag_spherlaguerre_sampling(nodes, weights, tau, P);

	int bandlimit_p = P;
	int bandlimit_l = L;

	offset_lmnp = 0;
	offset = 0;
	for (jp = parameters->J_min_p; jp <= J_p; jp++){
		if (!parameters->upsample) {
			bandlimit_p = MIN(flaglet_radial_bandlimit(jp, parameters), P);
			flag_spherlaguerre_sampling(nodes, weights, tau, bandlimit_p);
		}
		for (jl = parameters->J_min_l; jl <= J_l; jl++){
			if (!parameters->upsample) {
	            bandlimit_l = MIN(flaglet_angular_bandlimit(jl, parameters), L);
	            so3_parameters.L = bandlimit_l;
	            Nj = MIN(N,bandlimit_l);
	            Nj += (Nj+N)%2; // ensure N and Nj are both even or both odd
	            so3_parameters.N = Nj;
	        }
	        so3_parameters.L0 = ceil(pow(parameters->B_l, jl-1));
			lmn_size = so3_sampling_flmn_size(&so3_parameters);
			f_size = so3_sampling_f_size(&so3_parameters);
			f_wav_p = (complex double*)calloc( f_size * bandlimit_p, sizeof(complex double));
	        flag_spherlaguerre_mapped_synthesis_adjoint(f_wav_p, f_wav + offset, nodes, bandlimit_p, tau, bandlimit_p, f_size);
	        for (p = 0; p < bandlimit_p; p++){
		        so3_adjoint_inverse_direct(
		        	f_wav_lmnp + offset_lmnp,
		            f_wav_p + p * f_size,
		            &so3_parameters
		        );
		        offset_lmnp += lmn_size;
		    }
			offset += f_size * bandlimit_p;
	        free(f_wav_p);
		}
	}

	flag_adjoint_synthesis(f_scal_lmp, f_scal, nodes, P, L, tau, P, spin);
	//TODO make this function adjoint!!!
	flaglet_synthesis_lmnp(flmp, f_wav_lmnp, f_scal_lmp, wav_lmp, scal_lmp, parameters);

	flag_adjoint_analysis(f, flmp, L, tau, P, spin);

	free(nodes);
	free(weights);
	free(f_wav_lmnp);
	free(f_scal_lmp);
	free(flmp);
}




void flaglet_analysis_lmnp_adjoint(complex double *flmp, const complex double *f_wav_lmnp, const complex double *f_scal_lmp, const complex double *wav_lmp, const double *scal_lmp, const flaglet_parameters_t *parameters)
{
	int offset, jl, jp, l, m, n, p, lmn_size, lm_ind, ln_ind, indlmn, indjjlnp, indlmp;
	complex double psi;
	int L = parameters->L;
	int N = parameters->N;
	int P = parameters->P;
	int J_l = flaglet_j_max(L, parameters->B_l);
	int J_p = flaglet_j_max(P, parameters->B_p);
	so3_parameters_t so3_parameters = {};
	fill_so3_angular_parameters(&so3_parameters, parameters);
	int bandlimit_p = P;
	int bandlimit_l = L;
	int Nj = N;

	offset = 0;
	for (jp = parameters->J_min_p; jp <= J_p; jp++){
        if (!parameters->upsample)
        	bandlimit_p = MIN(flaglet_radial_bandlimit(jp, parameters), P);
		for (jl = parameters->J_min_l; jl <= J_l; jl++){
        	if (!parameters->upsample)
        	{
            	bandlimit_l = MIN(flaglet_angular_bandlimit(jl, parameters), L);
            	so3_parameters.L = bandlimit_l;
            	Nj = MIN(N,bandlimit_l);
            	Nj += (Nj+N)%2;
            	so3_parameters.N = Nj;
        	}
			lmn_size = so3_sampling_flmn_size(&so3_parameters);
			for (p = 0; p < bandlimit_p; p++){
	       		for (n = -Nj+1; n < Nj; n+=2){
					for (l = MAX(ABS(parameters->spin), ABS(n)); l < bandlimit_l; l++){
		       			ssht_sampling_elm2ind(&ln_ind, l, n);
		       			indjjlnp = jp * (J_l + 1) * L * L * P   +  jl * L * L * P + p * L * L + ln_ind;
						psi = 8*PI*PI/(2*l+1) * wav_lmp[indjjlnp];
						for (m = -l; m <= l ; m++){
	       					ssht_sampling_elm2ind(&lm_ind, l, m);
							indlmp = p * L * L + lm_ind;
							so3_sampling_elmn2ind(&indlmn, l, m, n, &so3_parameters);
							flmp[indlmp] += f_wav_lmnp[offset + p * lmn_size + indlmn] * psi ;
						}
					}
				}
			}
			offset += lmn_size * bandlimit_p;
		}
	}

	for (p = 0; p < P; p++){
		for (l = ABS(parameters->spin); l < L; l++){
			for (m = -l; m <= l ; m++){
				ssht_sampling_elm2ind(&lm_ind, l, m);
				indlmp = p * L * L  +  lm_ind;
				flmp[indlmp] += f_scal_lmp[indlmp] * sqrt((4.0*PI)/(2.0*l+1.0)) * scal_lmp[p*L+l] ;
			}
		}
	}
}

















