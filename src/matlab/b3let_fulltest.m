% b3let_fulltest - Run all tests
%
% B3LET package to perform Wavelets on the Solid Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICEPSE.txt for license details

clear all;
close all;

% Main parameters
R = 1.0;
L = 16;
P = 16;
B_l = 3;
B_p = 2;
J_min_l = 0;
J_min_p = 1;

J_l = ceil(log(L) ./ log(B_l))
J_p = ceil(log(P) ./ log(B_p))

% Checks tilling of harmonic space for axysimmetric wavelets
[kappa kappa0] = b3let_axisym_tilling(B_l, B_p, L, P, J_min_l, J_min_p);
error_on_axisym_tilling = b3let_check_axisym_tilling(kappa, kappa0, L, P, B_l, B_p)

% Generate random 3D FLAG decomposition
flmn = zeros(P, L^2);
flmn = rand(size(flmn)) + sqrt(-1)*rand(size(flmn));
flmn = 2.*(flmn - (1+sqrt(-1))./2);

% Generate the corresponding field
f = flag_synthesis(flmn);

% Test exactness of 3D Wavelet transform - multiresolution
[f_wav, f_scal] = b3let_axisym_analysis(f);
f_rec = b3let_axisym_synthesis(f_wav, f_scal);
default = max(max(max(abs(f-f_rec))))

% Test exactness of 3D Wavelet transform - full resolution
[f_wav, f_scal] = b3let_axisym_analysis(f, 'Downsample', false);
f_rec = b3let_axisym_synthesis(f_wav, f_scal, 'Downsample', false);
default_fullresolution = max(max(max(abs(f-f_rec))))

% Test exactness of 3D Wavelet transform - multiresolution
[f_wav, f_scal] = b3let_axisym_analysis(f, 'Downsample', true);
f_rec = b3let_axisym_synthesis(f_wav, f_scal, 'Downsample', true);
default_multiresolution = max(max(max(abs(f-f_rec))))

% Test exactness of 3D Wavelet transform - custom values
[f_wav, f_scal] = b3let_axisym_analysis(f, 'B_l', B_l, 'B_p', B_p, 'L', L, 'P', P, 'J_min_l', J_min_l, 'J_min_p', J_min_p, 'R', R);
f_rec = b3let_axisym_synthesis(f_wav, f_scal, 'B_l', B_l, 'B_p', B_p, 'L', L, 'P', P, 'J_min_l', J_min_l, 'J_min_p', J_min_p, 'R', R);
custom = max(max(max(abs(f-f_rec))))

% Impose reality on flms.
for en = 1:P
   for el = 0:L-1
      ind = el*el + el + 1;
      flmn(en,ind) = real(flmn(en,ind));
      for m = 1:el
         ind_pm = el*el + el + m + 1;
         ind_nm = el*el + el - m + 1;
         flmn(en,ind_nm) = (-1)^m * conj(flmn(en,ind_pm));
      end  
   end
end

% Generate the corresponding field
f_real = flag_synthesis(flmn, 'reality', true);

% Test exactness of 3D Wavelet transform - multiresolution
[f_wav_real, f_scal_real] = b3let_axisym_analysis(f_real, 'reality', true);
f_real_rec = b3let_axisym_synthesis(f_wav_real, f_scal_real, 'reality', true);
real_default = max(max(max(abs(f_real-f_real_rec))))

% Test exactness of 3D Wavelet transform - full resolution
[f_wav_real, f_scal_real] = b3let_axisym_analysis(f_real, 'Downsample', false, 'reality', true);
f_rec = b3let_axisym_synthesis(f_wav_real, f_scal_real, 'Downsample', false, 'reality', true);
real_default_fullresolution = max(max(max(abs(f_real-f_real_rec))))

% Test exactness of 3D Wavelet transform - multiresolution
[f_wav_real, f_scal_real] = b3let_axisym_analysis(f_real, 'Downsample', true, 'reality', true);
f_rec = b3let_axisym_synthesis(f_wav_real, f_scal_real, 'Downsample', true, 'reality', true);
real_default_multiresolution = max(max(max(abs(f_real-f_real_rec))))

% Test exactness of 3D Wavelet transform
[f_wav_real, f_scal_real] = b3let_axisym_analysis(f_real, 'B_l', B_l, 'B_p', B_p, 'L', L, 'P', P, 'J_min_l', J_min_l, 'J_min_p', J_min_p, 'R', R, 'reality', true);
f_real_rec = b3let_axisym_synthesis(f_wav_real, f_scal_real, 'B_l', B_l, 'B_p', B_p, 'L', L, 'P', P, 'J_min_l', J_min_l, 'J_min_p', J_min_p, 'R', R, 'reality', true);
real_custom = max(max(max(abs(f_real-f_real_rec))))
