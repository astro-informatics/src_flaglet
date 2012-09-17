function [f_wav, f_scal] = flaglet_axisym_analysis(f, varargin)

% flaglet_axisym_analysis 
% Compute axisymmetric wavelet transform, output in pixel space.
%
% Default usage :
%
%   [f_wav, f_scal] = flaglet_axisym_analysis(f, <options>)
%
% f is the input field -- MW sampling,
% f_wav contains the output wavelet contributions,
% f_scal contains the output scaling contributions,
% B_l is the wavelet parameter for angular space,
% B_p is the wavelet parameter for radial space,
% L is the angular band-limit,
% P is the radial band-limit,
% J_min_l the first angular wavelet scale to use,
% J_min_p the first radial wavelet scale to use,
% R is the radial boundary-limit.
%
% Options :
%  'B_l'               = { Dilation factor; B_l > 1 (default=2) }
%  'B_p'               = { Dilation factor; B_p > 1 (default=2) }
%  'L'               = { Angular harmonic band-limit; L > 1 (default=guessed) }
%  'P'               = { Radial harmonic band-limit; P > 1 (default=guessed) }
%  'J_min_l'           = { Minimum needlet scale to consider;
%                        0 <= J_min_l < log_B_l(L) (default=0) }
%  'J_min_p'           = { Minimum needlet scale to consider;
%                        0 <= J_min_p < log_B_p(P) (default=0) }
%  'R'               = { Radial boundary; R > 0 (default=1.0) }
%  'Downsample'      = { true        [multiresolution algorithm (default)],
%                        false       [full resolution wavelets] }
%  'Reality'         = { false        [do not assume f real (default)],
%                        true         [assume f real (improves performance)] }
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
p.addParamValue('L', Lguessed, @isnumeric);   
p.addParamValue('P', Pguessed, @isnumeric);    
p.addParamValue('J_min_l', 0, @isnumeric);   
p.addParamValue('J_min_p', 0, @isnumeric); 
p.addParamValue('R', 1.0, @isnumeric); 
p.addParamValue('Downsample', true, @islogical);
p.addParamValue('Reality', false, @islogical);
p.parse(f, varargin{:});
args = p.Results;

P = args.P;
L = args.L;
J_min_l = args.J_min_l;
J_min_p = args.J_min_p;
f_vec = zeros(P, L*(2*L-1));
for p = 1:args.P
    temp(:,:) = f(p,:,:);
    f_vec(p,:) = flag_mw_arr2vec( temp );
end

[f_wav_vec, f_scal_vec] = flaglet_axisym_analysis_mex(f_vec, args.B_l, args.B_p, args.L, args.P, args.J_min_l, args.J_min_p, args.R, args.Reality, args.Downsample);

J_l = ceil(log(L) ./ log(args.B_l));
J_p = ceil(log(P) ./ log(args.B_p));
f_wav = cell(J_l+1-J_min_l, J_p+1-J_min_p);

offset = 0;
for jp = J_min_p:J_p
    if args.Downsample == true
        band_limit_p = min([ ceil(args.B_p^(jp+1)) P ]);
     else
        band_limit_p = P;
    end
    for jl = J_min_l:J_l
        if args.Downsample == true
            band_limit_l = min([ ceil(args.B_l^(jl+1)) L ]);
        else
            band_limit_l = L;
        end
        temp = zeros(band_limit_p, band_limit_l, 2*band_limit_l-1);
        for r = 0:band_limit_p-1
            for t = 0:band_limit_l-1
                for p = 0:2*band_limit_l-2
                    ind = offset + r * band_limit_l *( 2 * band_limit_l - 1 ) + t * ( 2 * band_limit_l - 1) + p + 1;
                    temp(r+1,t+1,p+1) = f_wav_vec(1,ind);
                end
            end
        end
        f_wav{jl+1-J_min_l, jp+1-J_min_p} = temp;
        offset = offset + band_limit_l * (2*band_limit_l-1) * band_limit_p;
    end
end

f_scal = zeros(P, L, (2*L-1));
for p = 1:P
    temp = f_scal_vec(p,:);
    size(flag_mw_vec2arr( temp ));
    f_scal(p,:,:) = flag_mw_vec2arr( temp );
end

end