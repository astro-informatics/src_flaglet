function error_on_axisym_tiling = flaglet_check_axisym_tiling(kappa, kappa0, L, N, B_l, B_n)

% flaglet_check_axisym_tiling - Checks exactness of the tiling.
% -- Axisymmetric wavelets on the solid sphere.
%
% B3LET package to perform Wavelets transform on the Solid Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

J_l = ceil(log(L) ./ log(B_l));
J_n = ceil(log(N) ./ log(B_n));

identity = kappa0.^2;
for jl = 0:J_l
    for jn = 0:J_n
        temp = kappa{jl+1,jn+1};
        identity(:,:) = identity(:,:) + temp.^2;
    end
end

error_on_axisym_tiling = 0;
for l=1:L
    for n=1:N
        error_on_axisym_tiling = error_on_axisym_tiling + identity(n,l) - 1.0; 
    end
end

end