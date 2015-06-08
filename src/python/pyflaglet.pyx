
# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# if you want to use the Numpy-C-API from Cython
np.import_array()

#----------------------------------------------------------------------------------------------------#

cdef extern from "flaglet.h":

	ctypedef struct flaglet_parameters_t:
		int J_min_l
		int J_min_p
		int B_l
		int B_p
		int L
		int N
		int P
		int spin
		int upsample
		int reality
		double R

	void fill_so3_angular_parameters(so3_parameters_t *so3_parameters, const flaglet_parameters_t *parameters);

	int flaglet_radial_bandlimit(int jp, const flaglet_parameters_t *parameters);

	int flaglet_angular_bandlimit(int jl, const flaglet_parameters_t *parameters);

	void flaglet_analysis(double complex *f_wav, double complex *f_scal, const double complex *f, const flaglet_parameters_t *parameters);

	void flaglet_synthesis(double complex *f, const double complex *f_wav, const double complex *f_scal, const flaglet_parameters_t *parameters);

	void flaglet_allocate_f_wav(double complex **f_wav, double complex **f_scal, const flaglet_parameters_t *parameters);

#----------------------------------------------------------------------------------------------------#

cdef extern from "flag.h":

	void flag_core_analysis(double complex *flmn,
		const double complex *f,
		double R, int L, int N);

	void flag_core_synthesis(double complex *f,
		const double complex *flmn,
		const double *nodes,
		int Nnodes, int L, int N);

	void flag_sampling_allocate(double **rs, double **thetas, double **phis, double **laguweights, double R, int L, int N);

	void flag_sampling(double *rs, double *thetas, double *phis, double *laguweights, double R, int L, int N);

	void flag_spherlaguerre_sampling(double *nodes, double *weights, double R, int P)

	void flag_spherlaguerre_synthesis_gen(double *f, const double *fn, const double *nodes, int Nnodes, int N, int alpha)

#----------------------------------------------------------------------------------------------------#

cdef extern from "s2let.h":

	ctypedef enum s2let_sampling_t:
		S2LET_SAMPLING_MW
	ctypedef enum s2let_wav_norm_t:
		S2LET_WAV_NORM_DEFAULT, S2LET_WAV_NORM_SPIN_LOWERED

	ctypedef struct s2let_parameters_t:
		int J_min
		int B
		int L
		int N
		int upsample
		int spin
		ssht_dl_method_t dl_method
		s2let_wav_norm_t normalization;
		s2let_sampling_t sampling_scheme;
		int original_spin
		int reality
		int verbosity


#----------------------------------------------------------------------------------------------------#

cdef extern from "ssht.h":

	ctypedef enum ssht_dl_method_t:
		SSHT_DL_RISBO, SSHT_DL_TRAPANI

	int ssht_sampling_mw_ntheta(int L);
	int ssht_sampling_mw_nphi(int L);

#----------------------------------------------------------------------------------------------------#

cdef extern from "so3.h":

	ctypedef struct so3_parameters_t:
		int L0
		int L
		int N
		ssht_dl_method_t dl_method

	int s2let_j_max(s2let_parameters_t *parameters);

	int so3_sampling_f_size(const so3_parameters_t *parameters);

#----------------------------------------------------------------------------------------------------#

cdef extern from "stdlib.h":
	void free(void* ptr)

#----------------------------------------------------------------------------------------------------#

def j_max(B, L, J_min):
	cdef s2let_parameters_t parameters = {};
	parameters.B = B;
	parameters.L = L;
	parameters.J_min = J_min;
	return s2let_j_max(&parameters);

#----------------------------------------------------------------------------------------------------#

def sampling(R, L, N):

	ntheta = ssht_sampling_mw_ntheta(L);
	nphi = ssht_sampling_mw_nphi(L);
	thetas = np.empty([ntheta,], dtype=float)
	phis = np.empty([nphi,], dtype=float)
	rs = np.empty([N,], dtype=float)
	laguweights = np.empty([N,], dtype=float)

	flag_sampling(
		<double*> np.PyArray_DATA(rs),
		<double*> np.PyArray_DATA(thetas),
		<double*> np.PyArray_DATA(phis),
		<double*> np.PyArray_DATA(laguweights),
		R, L, N)

	return rs, thetas, phis

#----------------------------------------------------------------------------------------------------#

def spherlaguerre_synthesis(np.ndarray[double, ndim=1, mode="c"] fn not None,
	np.ndarray[double, ndim=1, mode="c"] nodes not None, int alpha):

	Nnodes = len(nodes)
	N = len(fn)
	f = np.zeros([Nnodes], dtype=float)

	flag_spherlaguerre_synthesis_gen(
		<double*> np.PyArray_DATA(f),
		<double*> np.PyArray_DATA(fn),
		<double*> np.PyArray_DATA(nodes),
		Nnodes, N, alpha);

	return f

#----------------------------------------------------------------------------------------------------#

def flag_analysis(np.ndarray[double complex, ndim=3, mode="c"] f not None,
	double R, int L, int N):

	flmn = np.zeros([N*L*L], dtype=complex)

	flag_core_analysis(
		<double complex*> np.PyArray_DATA(flmn),
		<double complex*> np.PyArray_DATA(f.ravel()),
		R, L, N);

	return flmn.reshape((N, L*L))

#----------------------------------------------------------------------------------------------------#

def flag_synthesis(np.ndarray[double complex, ndim=2, mode="c"] flmn not None,
	np.ndarray[double, ndim=1, mode="c"] nodes not None,
	int L, int N):

	f = np.zeros([N*L*(2*L-1)], dtype=complex)

	Nnodes = len(nodes)

	flag_core_synthesis(
		<double complex*> np.PyArray_DATA(f),
		<double complex*> np.PyArray_DATA(flmn.ravel()),
		<double*> np.PyArray_DATA(nodes),
		Nnodes, L, N);

	return f.reshape((N,L,2*L-1))

#----------------------------------------------------------------------------------------------------#

def pyflaglet_analysis(
	np.ndarray[double complex, ndim=3, mode="c"] f not None,
	R, L, N, P, J_min_l, J_min_p, B_l, B_p, spin, upsample):

	cdef flaglet_parameters_t parameters = {};
	parameters.N = N;
	parameters.P = P;
	parameters.R = R;
	parameters.L = L;
	parameters.J_min_l = J_min_l;
	parameters.J_min_p = J_min_p;
	parameters.B_l = B_l;
	parameters.B_p = B_p;
	parameters.spin = spin;
	parameters.upsample = upsample;

	J_l = j_max(B_l, L, J_min_l)
	J_p = j_max(B_p, P, J_min_p)

	offset, bandlimit_l, bandlimit_p, n_elem_wav = flaglet_wav_ind(N, J_l, J_p, L, N, P, B_l, B_p, J_min_l, J_min_p, upsample)
	total = offset + bandlimit_l*(2*bandlimit_l-1)*bandlimit_p
	f_wav = np.empty([total], dtype=complex)
	f_scal = np.empty([L*(2*L-1)*P], dtype=complex)

	flaglet_analysis(
		<double complex*> np.PyArray_DATA(f_wav),
		<double complex*> np.PyArray_DATA(f_scal),
		<double complex*> np.PyArray_DATA(f.ravel()),
		&parameters);


	return f_wav, f_scal

#----------------------------------------------------------------------------------------------------#

def pyflaglet_synthesis(
	np.ndarray[double complex, ndim=1, mode="c"] f_wav not None,
	np.ndarray[double complex, ndim=1, mode="c"] f_scal not None,
	R, L, N, P, J_min_l, J_min_p, B_l, B_p, spin, upsample):

	cdef flaglet_parameters_t parameters = {};
	parameters.N = N;
	parameters.P = P;
	parameters.L = L;
	parameters.R = R;
	parameters.J_min_l = J_min_l;
	parameters.J_min_p = J_min_p;
	parameters.B_l = B_l;
	parameters.B_p = B_p;
	parameters.spin = spin;
	parameters.upsample = upsample;

	f = np.empty([L*(2*L-1)*P], dtype=complex)

	flaglet_synthesis(
		<double complex*> np.PyArray_DATA(f),
		<double complex*> np.PyArray_DATA(f_wav),
		<double complex*> np.PyArray_DATA(f_scal),
		&parameters);

	return f.reshape((P,L,2*L-1))

#----------------------------------------------------------------------------------------------------#

def flaglet_wav_ind(n, jl, jp, L, N, P, B_l, B_p, J_min_l, J_min_p, upsample):

	cdef flaglet_parameters_t parameters = {};
	parameters.N = N;
	parameters.P = P;
	parameters.L = L;
	parameters.J_min_l = J_min_l;
	parameters.J_min_p = J_min_p;
	parameters.B_l = B_l;
	parameters.B_p = B_p;
	parameters.upsample = upsample;

	cdef so3_parameters_t so3_parameters = {};
	fill_so3_angular_parameters(&so3_parameters, &parameters);

	J_l = j_max(B_l, L, J_min_l)
	J_p = j_max(B_p, P, J_min_p)

	bandlimit_p = P
	bandlimit_l = L

	offset = 0
	for jprime_p from J_min_p <= jprime_p <= J_p:
		if upsample == 0:
			bandlimit_p = min(flaglet_radial_bandlimit(jprime_p, &parameters), P);
		for jprime_l from J_min_l <= jprime_l <= J_l:
			if upsample == 0:
				bandlimit_l = min(flaglet_angular_bandlimit(jprime_l, &parameters), L);
			so3_parameters.L = bandlimit_l;
			for en from 1 <= en <= N:
				offset += so3_sampling_f_size(&so3_parameters) * bandlimit_p;
				if en == n and jprime_l == jl and jprime_p == jp:
					n_elem_wav = so3_sampling_f_size(&so3_parameters) * bandlimit_p;
					offset -= so3_sampling_f_size(&so3_parameters) * bandlimit_p;
					return offset, bandlimit_l, bandlimit_p, n_elem_wav

	return 0

#----------------------------------------------------------------------------------------------------#
