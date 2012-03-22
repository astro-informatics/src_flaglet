function [f_wav, f_scal] = b3let_axisym_analysis(f, B_l, B_n, L, N, J_min_l, J_min_n, varargin)

% b3let_axisym_analysis 
% Compute axisymmetric wavelet transform, output in pixel space.
%
% Default usage :
%
%   [f_wav, f_scal] = b3let_axisym_analysis(f, B_l, B_n, L, N, J_min_l, J_min_n, <options>)
%
% f is the input field -- MW sampling,
% f_wav contains the output wavelet contributions,
% f_scal contains the output scaling contributions,
% B_l is the wavelet parameter for angular space,
% B_n is the wavelet parameter for radial space,
% L is the angular band-limit,
% N is the radial band-limit,
% J_min_l the first angular wavelet scale to use,
% J_min_n the first radial wavelet scale to use.
%
% Option :
%  'Reality'         = { false        [do not assume f real (default)],
%                        true         [assume f real (improves performance)] }
%
% B3LET package to perform Wavelets transform on the Solid Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

p = inputParser;
p.addRequired('f', @isnumeric); 
p.addRequired('B_l', @isnumeric);       
p.addRequired('B_n', @isnumeric);          
p.addRequired('L', @isnumeric);   
p.addRequired('N', @isnumeric);    
p.addRequired('J_min_l', @isnumeric);   
p.addRequired('J_min_n', @isnumeric); 
p.addParamValue('Reality', false, @islogical);
p.parse(f, B_l, B_n, L, N, J_min_l, J_min_n, varargin{:});
args = p.Results;

f_vec = zeros(N, L*(2*L-1));
for n = 1:N
    temp(:,:) = f(n,:,:);
    f_vec(n,:) = flag_mw_arr2vec( temp );
end

[f_wav_vec, f_scal_vec] = b3let_axisym_analysis_mex(f_vec, B_l, B_n, L, N, J_min_l, J_min_n, args.Reality);

J_l = ceil(log(L) ./ log(B_l));
J_n = ceil(log(N) ./ log(B_n));
f_wav = cell(J_l+1, J_n+1);
for jl = 0:J_l
    for jn = 0:J_n
        temp = zeros(N, L, 2*L-1);
        for n = 0:N-1
            for t = 0:L-1
                for p = 0:2*L-2
                    ind = jn*(J_l+1)*L*(2*L-1)*N + jl*L*(2*L-1)*N + n*L*(2*L-1) + t*(2*L-1) + p + 1;
                    temp(n+1,t+1,p+1) = f_wav_vec(1,ind);
                end
            end
        end
        f_wav{jl+1, jn+1} = temp;
    end
end

f_scal = zeros(N, L, (2*L-1));
for n = 1:N
    temp = f_scal_vec(n,:);
    size(temp)
    size(flag_mw_vec2arr( temp ));
    f_scal(n,:,:) = flag_mw_vec2arr( temp );
end

end