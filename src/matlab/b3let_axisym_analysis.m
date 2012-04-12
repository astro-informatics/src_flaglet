function [f_wav, f_scal] = b3let_axisym_analysis(f, varargin)

% b3let_axisym_analysis 
% Compute axisymmetric wavelet transform, output in pixel space.
%
% Default usage :
%
%   [f_wav, f_scal] = b3let_axisym_analysis(f, <options>)
%
% f is the input field -- MW sampling,
% f_wav contains the output wavelet contributions,
% f_scal contains the output scaling contributions,
% B_l is the wavelet parameter for angular space,
% B_n is the wavelet parameter for radial space,
% L is the angular band-limit,
% N is the radial band-limit,
% J_min_l the first angular wavelet scale to use,
% J_min_n the first radial wavelet scale to use,
% R is the radial boundary-limit.
%
% Options :
%  'B_l'               = { Dilation factor; B_l > 1 (default=2) }
%  'B_n'               = { Dilation factor; B_n > 1 (default=2) }
%  'L'               = { Angular harmonic band-limit; L > 1 (default=guessed) }
%  'N'               = { Radial harmonic band-limit; N > 1 (default=guessed) }
%  'J_min_l'           = { Minimum needlet scale to consider;
%                        0 <= J_min_l < log_B_l(L) (default=0) }
%  'J_min_n'           = { Minimum needlet scale to consider;
%                        0 <= J_min_n < log_B_n(N) (default=0) }
%  'R'               = { Radial boundary; R > 0 (default=1.0) }
%  'Reality'         = { false        [do not assume f real (default)],
%                        true         [assume f real (improves performance)] }
%
% B3LET package to perform Wavelets transform on the Solid Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

sz = size(f);
Nguessed = sz(1)-1;
Lguessed = sz(2);

p = inputParser;
p.addRequired('f', @isnumeric); 
p.addParamValue('B_l', 2, @isnumeric);       
p.addParamValue('B_n', 2, @isnumeric);          
p.addParamValue('L', Lguessed, @isnumeric);   
p.addParamValue('N', Nguessed, @isnumeric);    
p.addParamValue('J_min_l', 0, @isnumeric);   
p.addParamValue('J_min_n', 0, @isnumeric); 
p.addParamValue('R', 1.0, @isnumeric); 
p.addParamValue('Reality', false, @islogical);
p.parse(f, varargin{:});
args = p.Results;

N = args.N;
L = args.L;
f_vec = zeros(N+1, L*(2*L-1));
for n = 1:args.N+1
    temp(:,:) = f(n,:,:);
    f_vec(n,:) = flag_mw_arr2vec( temp );
end

[f_wav_vec, f_scal_vec] = b3let_axisym_analysis_mex(f_vec, args.B_l, args.B_n, args.L, args.N, args.J_min_l, args.J_min_n, args.R, args.Reality);

J_l = ceil(log(L) ./ log(args.B_l));
J_n = ceil(log(N) ./ log(args.B_n));
f_wav = cell(J_l+1, J_n+1);
for jl = 0:J_l
    for jn = 0:J_n
        temp = zeros(N+1, L, 2*L-1);
        for n = 0:N
            for t = 0:L-1
                for p = 0:2*L-2
                    ind = jn*(J_l+1)*L*(2*L-1)*(N+1) + jl*L*(2*L-1)*(N+1) + n*L*(2*L-1) + t*(2*L-1) + p + 1;
                    temp(n+1,t+1,p+1) = f_wav_vec(1,ind);
                end
            end
        end
        f_wav{jl+1, jn+1} = temp;
    end
end

f_scal = zeros(N+1, L, (2*L-1));
for n = 1:N+1
    temp = f_scal_vec(n,:);
    size(flag_mw_vec2arr( temp ));
    f_scal(n,:,:) = flag_mw_vec2arr( temp );
end

end