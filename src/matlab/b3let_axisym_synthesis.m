function f = b3let_axisym_synthesis(f_wav, f_scal, B_l, B_n, L, N, J_min_l, J_min_n, varargin)

% b3let_axisym_synthesis 
% Compute axisymmetric wavelet transform, output in pixel space.
%
% Default usage :
%
%   f = s2let_axisym_synthesis(f_wav, f_scal, B, L, J_min, <options>)
%
% f_wav contains the input wavelet contributions -- MW sampling,
% f_scal contains the input scaling contributions -- MW sampling,
% f is the output field -- MW sampling,
% B is the wavelet parameter,
% L is the angular band-limit,
% J_min the first wavelet to be used.
%
% Option :
%  'Reality'         = { false        [do not assume f real (default)],
%                        true         [assume f real (improves performance)] }
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

p = inputParser;
p.addRequired('f_wav'); 
p.addRequired('f_scal', @isnumeric); 
p.addRequired('B_l', @isnumeric);       
p.addRequired('B_n', @isnumeric);          
p.addRequired('L', @isnumeric);   
p.addRequired('N', @isnumeric);    
p.addRequired('J_min_l', @isnumeric);   
p.addRequired('J_min_n', @isnumeric); 
p.addParamValue('Reality', false, @islogical);
p.parse(f_wav, f_scal, B_l, B_n, L, N, J_min_l, J_min_n, varargin{:});
args = p.Results;

f_scal_vec = zeros(N, L*(2*L-1));
for n = 1:N
    temp(:,:) = f_scal(n,:,:);
    f_scal_vec(n,:) = flag_mw_arr2vec( temp );
end

J_l = ceil(log(L) ./ log(B_l));
J_n = ceil(log(N) ./ log(B_n));
f_wav_vec = zeros(1,(J_l+1)*(J_n+1)*L*(2*L-1)*N);
for jl = 0:J_l
    for jn = 0:J_n
        for n = 0:N-1
            for t = 0:L-1
                for p = 0:2*L-2
                    temp = f_wav{jl+1, jn+1};
                    ind = jn*(J_l+1)*L*(2*L-1)*N + jl*L*(2*L-1)*N + n*L*(2*L-1) + t*(2*L-1) + p + 1;
                    f_wav_vec(1,ind) = temp(n+1,t+1,p+1);
                end
            end
        end
    end
end

f_vec = b3let_axisym_synthesis_mex(f_wav_vec, f_scal_vec, B_l, B_n, L, N, J_min_l, J_min_n, args.Reality);
 
f = zeros(N, L, (2*L-1));
for n = 1:N
    temp = f_vec(n,:);
    f(n,:,:) = flag_mw_vec2arr( temp );
end

end