function [f_wav, f_scal] = flaglet_analysis(f, varargin)

% flaglet_analysis 
% Compute wavelet transform, output in pixel space.
%
% Default usage :
%
%   [f_wav, f_scal] = flaglet_analysis(f, <options>)
%
% f is the input field -- MW sampling,
% f_wav contains the output wavelet contributions,
% f_scal contains the output scaling contributions,
% B_l is the wavelet parameter for angular space,
% B_p is the wavelet parameter for radial space,
% L is the angular band-limit,
% P is the radial band-limit,
% N is the azimuthal band-limit,
% J_min_l the first angular wavelet scale to use,
% J_min_p the first radial wavelet scale to use,
% R is the radial boundary-limit.
%
% Options :
%  'B_l'               = { Dilation factor; B_l > 1 (default=2) }
%  'B_p'               = { Dilation factor; B_p > 1 (default=2) }
%  'N'               = { Azimuthal band-limit; N > 0 (default=1) }
%  'L'               = { Angular harmonic band-limit; L > 1 (default=guessed) }
%  'P'               = { Radial harmonic band-limit; P > 1 (default=guessed) }
%  'J_min_l'           = { Minimum needlet scale to consider;
%                        0 <= J_min_l < log_B_l(L) (default=0) }
%  'J_min_p'           = { Minimum needlet scale to consider;
%                        0 <= J_min_p < log_B_p(P) (default=0) }
%  'tau'               = { Radial boundary; R > 0 (default=1.0) }
%  'Downsample'      = { true        [multiresolution algorithm (default)],
%                        false       [full resolution wavelets] }
%
% flaglet package to perform Wavelets transform on the Solid Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

sz = size(f);
Pguessed = sz(1);
Lguessed = sz(2);

p = inputParser;
p.addRequired('f', @isnumeric); 
p.addParamValue('B_l', 2, @isnumeric);       
p.addParamValue('B_p', 2, @isnumeric);        
p.addParamValue('N', 1, @isnumeric);       
p.addParamValue('L', Lguessed, @isnumeric);   
p.addParamValue('P', Pguessed, @isnumeric);    
p.addParamValue('J_min_l', 0, @isnumeric);   
p.addParamValue('J_min_p', 0, @isnumeric); 
p.addParamValue('tau', 1.0, @isnumeric); 
p.addParamValue('Downsample', true, @islogical);
p.parse(f, varargin{:});
args = p.Results;

P = args.P;
L = args.L;
N = args.N;
J_min_l = args.J_min_l;
J_min_p = args.J_min_p;
f_vec = zeros(P, L*(2*L-1));
for p = 1:args.P
    temp(:,:) = f(p,:,:);
    f_vec(p,:) = flag_mw_arr2vec( temp );
end

[f_wav_vec, f_scal_vec] = flaglet_analysis_mex(f_vec, args.B_l, args.B_p, args.L, args.P, args.J_min_l, args.J_min_p, args.tau, args.N, args.Downsample);

J_l = ceil(log(L) ./ log(args.B_l));
J_p = ceil(log(P) ./ log(args.B_p));
f_wav = cell(J_l+1-J_min_l, J_p+1-J_min_p);

offset = 0;
for jp = J_min_p:J_p
    if args.Downsample == true
        band_limit_p = min([ s2let_bandlimit(jp,J_min_p,args.B_p,P) P ]);
     else
        band_limit_p = P;
    end
    for jl = J_min_l:J_l
        if args.Downsample == true
            band_limit_l = min([ s2let_bandlimit(jl,J_min_l,args.B_l,L) L ]);
        else
            band_limit_l = L;
        end
        temp = zeros(N, band_limit_p, band_limit_l, 2*band_limit_l-1);
        for r = 0:band_limit_p-1
            for n = 0:N-1
                for t = 0:band_limit_l-1
                    for p = 0:2*band_limit_l-2
                        ind = offset + r * N * band_limit_l *( 2 * band_limit_l - 1 ) + n * band_limit_l *( 2 * band_limit_l - 1 ) + t * ( 2 * band_limit_l - 1) + p + 1;
                        temp(n+1,r+1,t+1,p+1) = f_wav_vec(1,ind);
                    end
                end
            end
        end
        f_wav{jl+1-J_min_l, jp+1-J_min_p} = temp;
        offset = offset + band_limit_l * (2*band_limit_l-1) * band_limit_p * N;
    end
end

f_scal = zeros(P, L, (2*L-1));
for p = 1:P
    temp = f_scal_vec(p,:);
    size(flag_mw_vec2arr( temp ));
    f_scal(p,:,:) = flag_mw_vec2arr( temp );
end

end