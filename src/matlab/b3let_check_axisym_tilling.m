function error_on_axisym_tilling = b3let_check_axisym_tilling(kappa, kappa0, L, N, B_l, B_n)

% b3let_check_axisym_tilling - Checks exactness of the tilling.
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

error_on_axisym_tilling = 0;
for l=1:L
    for n=1:N
        error_on_axisym_tilling = error_on_axisym_tilling + identity(n,l) - 1.0; 
    end
end

end