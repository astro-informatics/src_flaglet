function [kappa kappa0] = b3let_axisym_tilling(B_l, B_n, L, N, J_min_l, J_min_n)

% b3let_axisym_tilling - Compute tilling in l-n harmonic space.
% -- Axisymmetric wavelets on the solid sphere.
%
% Default usage :
%
%   [kappa kappa0] = b3let_axisym_tilling(B_l, B_n, L, N, J_min_l, J_min_n)
%
% kappa is an array containing wavelet tilling functions,
% kappa0 contains the scaling function,
% B_l is the wavelet parameter for angular space,
% B_n is the wavelet parameter for radial space,
% L is the angular band-limit,
% N is the radial band-limit,
% J_min_l the first angular wavelet scale to use,
% J_min_n the first radial wavelet scale to use.
%
% B3LET package to perform Wavelets transform on the solid Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

p = inputParser;
p.addRequired('B_l', @isnumeric);       
p.addRequired('B_n', @isnumeric);          
p.addRequired('L', @isnumeric);   
p.addRequired('N', @isnumeric);    
p.addRequired('J_min_l', @isnumeric);   
p.addRequired('J_min_n', @isnumeric); 
p.parse(B_l, B_n, L, N, J_min_l, J_min_n);
args = p.Results;

[kappa_vec kappa0] = b3let_axisym_tilling_mex(B_l, B_n, L, N, J_min_l, J_min_n);

J_l = ceil(log(L) ./ log(B_l));
J_n = ceil(log(N) ./ log(B_n));

kappa = cell(J_l + 1, J_n + 1);
for jl = 0:J_l
    for jn = 0:J_n
        temp = zeros(N,L);
        for l = 0:L-1
            for n = 0:N-1
                temp(n+1,l+1) = kappa_vec(1,jn*(J_l+1)*L*N  + jl*L*N + n*L + l + 1);
            end
        end
        kappa{jl+1,jn+1} = temp;
    end   
end

end