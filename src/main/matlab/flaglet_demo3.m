% flaglet_demo3
%
% Plot wavelet decomposition of LSS dataset
%
% flaglet package to perform Wavelets transform on the Solid Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

L = 192 ;
P = 192 ;
Lplot = 192;
Pplot = 192;
minbound = 0;
maxbound = 30;
layer = 165;
zoomfactor = 1.9;

B_l = 7;
B_p = 6;
J_min_l = 1;
J_min_p = 2;
downsample = false;

J_l = min([J_min_l+1, ceil(log(L)/log(B_l))])
J_p = min([J_min_p+1, ceil(log(P)/log(B_p))])

% Load LSS data
R = 420;
load('horizon_01_30k')
f = flag_pixelizedeltas(horizon_01_30k, L, P, R);

%load('lss_flmp_312')
%f_lmp = zeros(P, L^2);
%f_lmp(:,:) = f_lmp_full(1:P, 1:L^2);
%disp('Flag_synthesis...')
%f = flag_synthesis(f_lmp, 'R', R, 'reality', true);

% subplot(1,2,1)
% flaglet_plot_f( f_ini )
% caxis([minbound maxbound])
% subplot(1,2,2)
% flaglet_plot_f( f )
% caxis([minbound maxbound])
 
disp('Wavelet analysis...')
[f_wav, f_scal] = flaglet_axisym_analysis(f, 'B_l', B_l, 'B_p', B_p, 'reality', true, 'downsample', downsample, 'J_min_l', J_min_l, 'J_min_p', J_min_p);

% nodes = slag_sampling(256, R);
% f_lmp_zeropadded = zeros(P, 256^2);
% f_lmp_zeropadded(1:P, 1:L^2) = f_lmp(:,:);
% disp('flag synthesis')
% f_oversampled = flag_synthesis(f_lmp_zeropadded, 'Nodes', nodes, 'reality', true);    


disp('Plot...')
figure('Position',[100 100 900 900]) 


J_l_range = J_min_l:J_l%[ 1 2 ];
J_n_range = J_min_p:J_p%[ 1 2 ];
nx = 3;%length(J_l)+1;
ny = 2;%length(J_n);

ind = 0;


ind = ind + 1;
subplot(ny, nx, ind);
flaglet_plot_f( f, 'L', Lplot, 'P', Pplot, 'layer', layer )
zoom(zoomfactor)
colormap(flipud(gray))
caxis([minbound maxbound])
title('Data')

ind = ind + 1;
subplot(ny, nx, ind);
flaglet_plot_f( f_scal, 'L', Lplot, 'P', Pplot, 'layer', layer )
zoom(zoomfactor)
colormap(flipud(gray))
caxis([minbound maxbound])
title('Scaling part')

for jl = J_l_range - J_min_l + 1
    for jn = J_n_range - J_min_p + 1
        ind = ind + 1; 
        subplot(ny, nx, ind);
        flaglet_plot_f( f_wav{jl, jn}, 'L', Lplot, 'P', Pplot, 'layer', layer )
        zoom(zoomfactor)
        colormap(flipud(gray))
        %title(['Wavelet : jl=',int2str(jl),' jn=',int2str(jn)])
        caxis([minbound maxbound])
    end
end

