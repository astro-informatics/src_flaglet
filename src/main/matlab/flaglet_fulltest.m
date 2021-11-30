% flaglet_fulltest - tauun all tests
%
% flaglet package to perform Wavelets on the Solid Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICEPSE.txt for license details

clear all;
close all;

% Main parameters
tau = 2.0;
N = 5;
L = 16;
P = 16;
B_l = 3;
B_p = 2;
J_min_l = 0;
J_min_p = 1;

J_l = ceil(log(L) ./ log(B_l))
J_p = ceil(log(P) ./ log(B_p))

% Checks tiling of harmonic space for axysimmetric wavelets
%[kappa kappa0] = flaglet_axisym_tiling(B_l, B_p, L, P, J_min_l, J_min_p);
%error_on_axisym_tiling = flaglet_check_axisym_tiling(kappa, kappa0, L, P, B_l, B_p)

% Generate random 3D FLAG decomposition
flmn = zeros(P, L^2);
flmn = rand(size(flmn)) + sqrt(-1)*rand(size(flmn));
flmn = 2.*(flmn - (1+sqrt(-1))./2);

% Generate the corresponding field
f = flag_synthesis(flmn);

f = flag_synthesis(flmn);

% Test exactness of 3D Wavelet transform - full resolution
[f_wav, f_scal] = flaglet_analysis(f, 'Downsample', false);
f_rec = flaglet_synthesis(f_wav, f_scal, 'Downsample', false);
default_fullresolution_axisym = max(max(max(abs(f-f_rec))))

% Test exactness of 3D Wavelet transform - multiresolution
[f_wav, f_scal] = flaglet_analysis(f, 'Downsample', true);
f_rec = flaglet_synthesis(f_wav, f_scal, 'Downsample', true);
default_multiresolution_axisym = max(max(max(abs(f-f_rec))))

% Test exactness of 3D Wavelet transform - full resolution
[f_wav, f_scal] = flaglet_analysis(f, 'Downsample', false, 'N', N);
f_rec = flaglet_synthesis(f_wav, f_scal, 'Downsample', false, 'N', N);
default_fullresolution_directional = max(max(max(abs(f-f_rec))))

% Test exactness of 3D Wavelet transform - multiresolution
[f_wav, f_scal] = flaglet_analysis(f, 'Downsample', true, 'N', N);
f_rec = flaglet_synthesis(f_wav, f_scal, 'Downsample', true, 'N', N);
default_multiresolution_directional = max(max(max(abs(f-f_rec))))


% Test exactness of 3D Wavelet transform - custom values
[f_wav, f_scal] = flaglet_analysis(f, 'B_l', B_l, 'B_p', B_p, 'L', L, 'P', P, 'J_min_l', J_min_l, 'J_min_p', J_min_p, 'tau', tau, 'N', N);
f_rec = flaglet_synthesis(f_wav, f_scal, 'B_l', B_l, 'B_p', B_p, 'L', L, 'P', P, 'J_min_l', J_min_l, 'J_min_p', J_min_p, 'tau', tau, 'N', N);
custom = max(max(max(abs(f-f_rec))))