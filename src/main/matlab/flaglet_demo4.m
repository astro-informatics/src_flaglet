% flaglet_demo4
%
% Plot wavelet decomposition of geophysics data
%
% B3LET package to perform Wavelets transform on the Solid Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

L = 50 ;
P = 50 ;
R = 3000;

B_l = 4;
B_p = 4;
J_min_l = 2;
J_min_p = 2;
downsample = true;

J_l = ceil(log(L)/log(B_l))
J_p = ceil(log(P)/log(B_p))

f = gen_geomodel(L, P);

disp('Flag_analysis...')
f_lmp = flag_analysis(f, 'R', R, 'reality', true);

% subplot(1,2,1)
% flaglet_plot_f( f_ini )
% caxis([minbound maxbound])
% subplot(1,2,2)
% flaglet_plot_f( f )
% caxis([minbound maxbound])
 
disp('Wavelet analysis...')
[f_wav, f_scal] = flaglet_axisym_analysis(f, 'B_l', B_l, 'B_p', B_p, 'reality', true, 'downsample', downsample, 'J_min_l', J_min_l, 'J_min_p', J_min_p);



disp('Plot...')
figure('Position',[100 100 600 900]) 

% Layout parameters
minbound = -1.0; 
maxbound = 1.0;
Lplot = 192;
Pplot = 192;
layer = 173;
zoomfactor = 2.2;

J_l = [ 1 2 ];
J_n = [ 1 2 ];
nx = 2;
ny = 3;


ind = 1;
subplot(ny, nx, ind);
flaglet_plot_f( f, 'L', Lplot, 'P', Pplot, 'layer', layer )
zoom(zoomfactor)
colormap(flipud(gray))
caxis([-1.4 1.4])
title('Data')

ind = 2;
subplot(ny, nx, ind);
flaglet_plot_f( f_scal, 'L', Lplot, 'P', Pplot, 'layer', layer )
zoom(zoomfactor)
colormap(flipud(gray))
caxis([-1.2 1.2])
title('Scaling part')

for jl = J_l
    for jn = J_n
        ind = ind + 1; 
        subplot(ny, nx, ind);
        flaglet_plot_f( f_wav{jl, jn}, 'L', Lplot, 'P', Pplot, 'layer', layer )
        zoom(zoomfactor)
        %colormap(flipud(gray))
        title(['Wavelet : jl=',int2str(jl),' jn=',int2str(jn)])
        caxis([minbound maxbound])
    end
end