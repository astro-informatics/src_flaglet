% flaglet_demo3
%
% Plot wavelet decomposition of LSS dataset
%
% flaglet package to perform Wavelets transform on the Solid Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

L = 128 ;
P = 128 ;
Lplot = 256;
Pplot = 256;

B_l = 4;
B_p = 4;
J_min_l = 2;
J_min_p = 3;
downsample = false;

J_l = ceil(log(L)/log(B_l))
J_p = ceil(log(P)/log(B_p))

% Load LSS data
load('lss_flmp_312')
R = 420;

f_lmp = zeros(P, L^2);
f_lmp(:,:) = f_lmp_full(1:P, 1:L^2);

disp('Flag_synthesis...')
f = flag_synthesis(f_lmp, 'R', R, 'reality', true);

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
figure('Position',[100 100 600 900]) 

minbound = 12;
maxbound = 100;
layer = 230;
zoomfactor = 1.9;

J_l = J_min_l:(J_l-1)%[ 1 2 ];
J_n = J_min_p:J_p%[ 1 2 ];
nx = 3%length(J_l)+1;
ny = 2%length(J_n);

ind = 1;
subplot(ny, nx, ind);
flaglet_plot_f( f, 'L', Lplot, 'P', Pplot, 'layer', layer )
zoom(zoomfactor)
colormap(flipud(gray))
caxis([minbound maxbound])
%title('Data')

ind = 2;
subplot(ny, nx, ind);
flaglet_plot_f( f_scal, 'L', Lplot, 'P', Pplot, 'layer', layer )
zoom(zoomfactor)
colormap(flipud(gray))
caxis([minbound maxbound])
%title('Scaling part')

for jl = J_l - J_min_l + 1
    for jn = J_n - J_min_p + 1
        ind = ind + 1; 
        subplot(ny, nx, ind);
        flaglet_plot_f( f_wav{jl, jn}, 'L', Lplot, 'P', Pplot, 'layer', layer )
        zoom(zoomfactor)
        colormap(flipud(gray))
        %title(['Wavelet : jl=',int2str(jl),' jn=',int2str(jn)])
        caxis([minbound maxbound])
    end
end

