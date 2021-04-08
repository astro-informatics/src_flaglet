import numpy as np 
from pytest import approx 

import pyflaglet as flaglet

def test_directional_transform():

	parameters = flaglet.flaglet_parameters(
		   B_l = 2, 
		     L = 32, 
	   J_min_l = 0, 
		     N = 4, 
		   B_p = 2, 
		     P = 32, 
	   J_min_p = 0, 
		  spin = 0, 
	  upsample = 1, 
	   reality = 0, 
		   tau = 1.0)

	f_wav_size = flaglet.flaglet_wav_dim(parameters)
	f_scal_size = flaglet.flaglet_scal_dim(parameters)
	f_size = flaglet.flaglet_f_dim(parameters)

	rng = np.random.default_rng()
	f = rng.normal(size=(f_size)) + 1j*rng.normal(size=(f_size))

	f_wav_input, f_scal_input = flaglet.flaglet_ana(f, parameters)
	f_input = flaglet.flaglet_syn(f_wav_input, f_scal_input, parameters)
	f_wav_output, f_scal_output = flaglet.flaglet_ana(f_input, parameters)
	f_output = flaglet.flaglet_syn(f_wav_output, f_scal_output, parameters)

	assert f_output == approx(f_input)
	assert f_wav_output == approx(f_wav_input)
	assert f_scal_output == approx(f_scal_input)

def test_directional_analysis_adjoint_transform():

	parameters = flaglet.flaglet_parameters(
		   B_l = 2, 
		     L = 32, 
	   J_min_l = 0, 
		     N = 4, 
		   B_p = 2, 
		     P = 32, 
	   J_min_p = 0, 
		  spin = 0, 
	  upsample = 1, 
	   reality = 0, 
		   tau = 1.0)

	f_wav_size = flaglet.flaglet_wav_dim(parameters)
	f_scal_size = flaglet.flaglet_scal_dim(parameters)
	f_size = flaglet.flaglet_f_dim(parameters)

	rng = np.random.default_rng()
	f_2 = rng.normal(size=(f_size)) + 1j*rng.normal(size=(f_size))
	f_1 = rng.normal(size=(f_size)) + 1j*rng.normal(size=(f_size))

	w1, s1 = flaglet.flaglet_ana(f_1, parameters)
	f_1 = flaglet.flaglet_syn(w1, s1, parameters)
	w2, s2 = flaglet.flaglet_ana(f_2, parameters)

	w1, s1 = flaglet.flaglet_ana(f_1, parameters)
	f_2 = flaglet.flaglet_ana_adjoint(w2, s2, parameters)
	
	assert np.abs(np.vdot(f_2, f_1)) == approx(np.abs(np.vdot(w1, w2) + np.vdot(s1, s2)))

def test_directional_synthesis_adjoint_transform():

	parameters = flaglet.flaglet_parameters(
		   B_l = 2, 
		     L = 32, 
	   J_min_l = 0, 
		     N = 4, 
		   B_p = 2, 
		     P = 32, 
	   J_min_p = 0, 
		  spin = 0, 
	  upsample = 1, 
	   reality = 0, 
		   tau = 1.0)

	f_wav_size = flaglet.flaglet_wav_dim(parameters)
	f_scal_size = flaglet.flaglet_scal_dim(parameters)
	f_size = flaglet.flaglet_f_dim(parameters)

	rng = np.random.default_rng()
	f_2 = rng.normal(size=(f_size)) + 1j*rng.normal(size=(f_size))
	f_1 = rng.normal(size=(f_size)) + 1j*rng.normal(size=(f_size))

	w1, s1 = flaglet.flaglet_ana(f_1, parameters)
	f_1 = flaglet.flaglet_syn(w1, s1, parameters)
	w2, s2 = flaglet.flaglet_ana(f_2, parameters)

	w1, s1 = flaglet.flaglet_syn_adjoint(f_1, parameters)
	f_2 = flaglet.flaglet_syn(w2, s2, parameters)
	
	assert np.abs(np.vdot(f_1, f_2)) == approx(np.abs(np.vdot(w2, w1) + np.vdot(s2, s1)))


