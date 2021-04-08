# cython: language_level=3

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# if you want to use the Numpy-C-API from Cython
np.import_array()

#----------------------------------------------------------------------------------------------------#


cdef extern from "flaglet/flaglet.h":

	# Parameter structure for flaglet
	ctypedef struct flaglet_parameters_t:
			int B_l	
			int L	
			int J_min_l	
			int N	
			int B_p	
			int P	
			int J_min_p	
			int spin	
			int upsample	
			int reality	
			double tau	

	# Dimensionality helper functions
	int flaglet_f_size(const flaglet_parameters_t *parameters);
	int flaglet_scal_size(const flaglet_parameters_t *parameters);
	int flaglet_wav_size(const flaglet_parameters_t *parameters);
	# Directional flaglet core functions
	void flaglet_analysis(double complex *f_wav, double complex *f_scal, const double complex *f, const flaglet_parameters_t *parameters);
	void flaglet_analysis_adjoint(double complex *f, const double complex *f_wav, const double complex *f_scal, const flaglet_parameters_t *parameters);
	void flaglet_synthesis(double complex *f, const double complex *f_wav, const double complex *f_scal, const flaglet_parameters_t *parameters);
	void flaglet_synthesis_adjoint(double complex *f_wav, double complex *f_scal, const double complex *f, const flaglet_parameters_t *parameters);

def flaglet_parameters(int B_l, int L, int J_min_l, int N, int B_p, int P, int J_min_p, int spin, int upsample, int reality, double tau):
	d = dict()
	d["B_l"] = B_l
	d["L"] = L
	d["J_min_l"] = J_min_l
	d["N"] = N
	d["B_p"] = B_p
	d["P"] = P
	d["J_min_p"] = J_min_p
	d["spin"] = spin
	d["upsample"] = upsample
	d["reality"] = reality
	d["tau"] = tau
	return d

def flaglet_wav_dim(dict parameters):
	cdef flaglet_parameters_t parameters_t = {}
	parameters_t.B_l = parameters["B_l"]
	parameters_t.L = parameters["L"]
	parameters_t.J_min_l = parameters["J_min_l"]
	parameters_t.N = parameters["N"]
	parameters_t.B_p = parameters["B_p"]
	parameters_t.P = parameters["P"]
	parameters_t.J_min_p = parameters["J_min_p"]
	parameters_t.spin = parameters["spin"]
	parameters_t.upsample = parameters["upsample"]
	parameters_t.reality = parameters["reality"]
	parameters_t.tau = parameters["tau"]
	return flaglet_wav_size(&parameters_t)

def flaglet_scal_dim(dict parameters):
	cdef flaglet_parameters_t parameters_t = {}
	parameters_t.B_l = parameters["B_l"]
	parameters_t.L = parameters["L"]
	parameters_t.J_min_l = parameters["J_min_l"]
	parameters_t.N = parameters["N"]
	parameters_t.B_p = parameters["B_p"]
	parameters_t.P = parameters["P"]
	parameters_t.J_min_p = parameters["J_min_p"]
	parameters_t.spin = parameters["spin"]
	parameters_t.upsample = parameters["upsample"]
	parameters_t.reality = parameters["reality"]
	parameters_t.tau = parameters["tau"]
	return flaglet_scal_size(&parameters_t)

def flaglet_f_dim(dict parameters):
	cdef flaglet_parameters_t parameters_t = {}
	parameters_t.B_l = parameters["B_l"]
	parameters_t.L = parameters["L"]
	parameters_t.J_min_l = parameters["J_min_l"]
	parameters_t.N = parameters["N"]
	parameters_t.B_p = parameters["B_p"]
	parameters_t.P = parameters["P"]
	parameters_t.J_min_p = parameters["J_min_p"]
	parameters_t.spin = parameters["spin"]
	parameters_t.upsample = parameters["upsample"]
	parameters_t.reality = parameters["reality"]
	parameters_t.tau = parameters["tau"]
	return flaglet_f_size(&parameters_t)


def flaglet_ana(np.ndarray[double complex, ndim=1, mode="c"] f not None, dict parameters):
	cdef flaglet_parameters_t parameters_t = {}
	parameters_t.B_l = parameters["B_l"]
	parameters_t.L = parameters["L"]
	parameters_t.J_min_l = parameters["J_min_l"]
	parameters_t.N = parameters["N"]
	parameters_t.B_p = parameters["B_p"]
	parameters_t.P = parameters["P"]
	parameters_t.J_min_p = parameters["J_min_p"]
	parameters_t.spin = parameters["spin"]
	parameters_t.upsample = parameters["upsample"]
	parameters_t.reality = parameters["reality"]
	parameters_t.tau = parameters["tau"]

	cdef int f_wav_size, f_scal_size
	f_wav_size = flaglet_wav_size(&parameters_t)
	f_scal_size = flaglet_scal_size(&parameters_t)

	f_wav = np.zeros([f_wav_size,], dtype=complex)
	f_scal = np.zeros([f_scal_size,], dtype=complex)

	flaglet_analysis(<double complex*> np.PyArray_DATA(f_wav), 
					 <double complex*> np.PyArray_DATA(f_scal),
					 <const double complex*> np.PyArray_DATA(f),
					 &parameters_t)
	return f_wav, f_scal

def flaglet_ana_adjoint(np.ndarray[double complex, ndim=1, mode="c"] f_wav not None, 
						np.ndarray[double complex, ndim=1, mode="c"] f_scal not None, 
						dict parameters):
	cdef flaglet_parameters_t parameters_t = {}
	parameters_t.B_l = parameters["B_l"]
	parameters_t.L = parameters["L"]
	parameters_t.J_min_l = parameters["J_min_l"]
	parameters_t.N = parameters["N"]
	parameters_t.B_p = parameters["B_p"]
	parameters_t.P = parameters["P"]
	parameters_t.J_min_p = parameters["J_min_p"]
	parameters_t.spin = parameters["spin"]
	parameters_t.upsample = parameters["upsample"]
	parameters_t.reality = parameters["reality"]
	parameters_t.tau = parameters["tau"]

	cdef int f_size = flaglet_f_size(&parameters_t)
	f = np.zeros([f_size,], dtype=complex)

	flaglet_analysis_adjoint(<double complex*> np.PyArray_DATA(f), 
					  		 <const double complex*> np.PyArray_DATA(f_wav),
					  		 <const double complex*> np.PyArray_DATA(f_scal),
					  		 &parameters_t)
	return f


def flaglet_syn(np.ndarray[double complex, ndim=1, mode="c"] f_wav not None, 
				np.ndarray[double complex, ndim=1, mode="c"] f_scal not None, 
				dict parameters):
	cdef flaglet_parameters_t parameters_t = {}
	parameters_t.B_l = parameters["B_l"]
	parameters_t.L = parameters["L"]
	parameters_t.J_min_l = parameters["J_min_l"]
	parameters_t.N = parameters["N"]
	parameters_t.B_p = parameters["B_p"]
	parameters_t.P = parameters["P"]
	parameters_t.J_min_p = parameters["J_min_p"]
	parameters_t.spin = parameters["spin"]
	parameters_t.upsample = parameters["upsample"]
	parameters_t.reality = parameters["reality"]
	parameters_t.tau = parameters["tau"]

	cdef int f_size = flaglet_f_size(&parameters_t)
	f = np.zeros([f_size,], dtype=complex)

	flaglet_synthesis(<double complex*> np.PyArray_DATA(f), 
					  <const double complex*> np.PyArray_DATA(f_wav),
					  <const double complex*> np.PyArray_DATA(f_scal),
					  &parameters_t)
	return f

def flaglet_syn_adjoint(np.ndarray[double complex, ndim=1, mode="c"] f not None, dict parameters):
	cdef flaglet_parameters_t parameters_t = {}
	parameters_t.B_l = parameters["B_l"]
	parameters_t.L = parameters["L"]
	parameters_t.J_min_l = parameters["J_min_l"]
	parameters_t.N = parameters["N"]
	parameters_t.B_p = parameters["B_p"]
	parameters_t.P = parameters["P"]
	parameters_t.J_min_p = parameters["J_min_p"]
	parameters_t.spin = parameters["spin"]
	parameters_t.upsample = parameters["upsample"]
	parameters_t.reality = parameters["reality"]
	parameters_t.tau = parameters["tau"]

	cdef int f_wav_size, f_scal_size
	f_wav_size = flaglet_wav_size(&parameters_t)
	f_scal_size = flaglet_scal_size(&parameters_t)

	f_wav = np.zeros([f_wav_size,], dtype=complex)
	f_scal = np.zeros([f_scal_size,], dtype=complex)

	flaglet_synthesis_adjoint(<double complex*> np.PyArray_DATA(f_wav), 
					 		  <double complex*> np.PyArray_DATA(f_scal),
					 		  <const double complex*> np.PyArray_DATA(f),
					 		  &parameters_t)
	return f_wav, f_scal






