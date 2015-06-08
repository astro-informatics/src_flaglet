function [kappa kappa0] = flaglet_tiling(B_l, B_p, L, P, N, spin, J_min_l, J_min_p)

% flaglet_tiling - Compute tiling in l-n harmonic space.
% -- Axisymmetric wavelets on the solid sphere.
%
% Default usage :
%
%   [kappa kappa0] = flaglet_axisym_tiling(B_l, B_n, L, N, J_min_l, J_min_n)
%
% kappa is an array containing wavelet tiling functions,
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
p.addRequired('B_p', @isnumeric);          
p.addRequired('L', @isnumeric);   
p.addRequired('P', @isnumeric);        
p.addRequired('N', @isnumeric);   
p.addRequired('spin', @isnumeric);    
p.addRequired('J_min_l', @isnumeric);   
p.addRequired('J_min_p', @isnumeric); 
p.parse(B_l, B_p, L, P, N, spin, J_min_l, J_min_p);
args = p.Results;

[kappa_vec kappa0_vec] = flaglet_tiling_mex(B_l, B_p, L, P, N, spin, J_min_l, J_min_p);

J_l = ceil(log(L) ./ log(B_l));
J_p = ceil(log(P) ./ log(B_p));

kappa0 = zeros(P,L*L);
for p = 0:P-1
    for l = 0:L-1
        for m = -l:l
            kappa0(p+1,l*l+l+m+1) = kappa0_vec(1,p*L*L +l*l+l+m+1);
        end
    end
end  

kappa = cell(J_l + 1, J_p + 1);
for jl = 0:J_l
    for jp = 0:J_p
        temp = zeros(P,L*L);
        for p = 0:P-1
            for l = 0:L-1
                for m = -l:l
                    temp(p+1,l*l+l+m+1) = kappa_vec(1, jp*(J_l+1)*L*L*P  + jl*L*L*P + p*L*L + l*l+l+m+1);
                end
            end
        end
        kappa{jl+1,jp+1} = temp;
    end   
end

end