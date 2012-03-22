% b3let_fulltest - Run all tests
%
% B3LET package to perform Wavelets on the Solid Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

clear all;
close all;

% Main parameters
L = 16
N = 16
B_l = 2
B_n = 2
J_min_l = 1
J_min_n = 1

J_l = ceil(log(L) ./ log(B_l))
J_n = ceil(log(N) ./ log(B_n))

% Checks tilling of harmonic space for axysimmetric wavelets
[kappa kappa0] = b3let_axisym_tilling(B_l, B_n, L, N, J_min_l, J_min_n);
error_on_axisym_tilling = b3let_check_axisym_tilling(kappa, kappa0, L, N, B_l, B_n)

% Generate random 3D FLAG decomposition
flmn = zeros(N, L^2);
flmn = rand(size(flmn)) + sqrt(-1)*rand(size(flmn));
flmn = 2.*(flmn - (1+sqrt(-1))./2);

% Generate the corresponding field
f = flag_synthesis(flmn, L, N);
% Test exactness of 3D Wavelet transform
[f_wav, f_scal] = b3let_axisym_analysis(f, B_l, B_n, L, N, J_min_l, J_min_n);

% Impose reality on flms.
for en = 1:N
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
f_real = flag_synthesis(flmn, L, N, 'reality', true);