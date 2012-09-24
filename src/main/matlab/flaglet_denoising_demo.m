% flaglet_denoising_demo
% Denoising example (geophysics data)

% Resolution parameters
L = 128 ;
P = 128 ;
R = 3000;

% Layout parameters
minbound = -1.4; 
maxbound = 1.4;
zoomfactor = 2.2;

% Generates geophysics data on the ball
f_ini = gen_geomodel(L, P);

disp('Analysis the signal (flag_analysis)...')
f_lmp = flag_analysis(f_ini, 'R', R, 'reality', true);
disp('Reconstruct the field (flag_synthesis)...')
f = flag_synthesis(f_lmp, 'R', R, 'reality', true);

% r = load('horizon_01_30k.mat');
% data = r.horizon_01_30k;
% 
% if max(data(:,2)) > max(data(:,3)) 
%     temp = data;
%     data(:,3) = temp(:,2);
%     data(:,2) = temp(:,3);
% end
% 
% R = 420;
% 
% disp('Pixelise data...')
% f_ini = pixelizedeltas( data, L, P, R );
% disp('Flag_analysis...')
% f_lmp = flag_analysis(f_ini, 'R', R, 'reality', true);
% disp('Flag_synthesis...')
% f = flag_synthesis(f_lmp, 'R', R, 'reality', true);

% Wavelet decomposition parameters
Downsample = true;
B_l = 3 ;
B_p = 3 ;
J_min_l = 0 ;
J_min_p = 0 ;
J_l = ceil(log(L) ./ log(B_l));
J_p = ceil(log(P) ./ log(B_p));

disp('Compute tilling of harmonic space...')
[kappa kappa0] = flaglet_axisym_tilling(B_l, B_p, L, P, J_min_l, J_min_p);

disp('Find noise covariance from SNR')
SNR_ini = 8.0;
f_power = sum(sum(abs(f_lmp).^2));
noise_power = f_power * ( 10.0^(-SNR_ini/10.0) );
factor = sum( ((1:P)./P).^2 ) * L * L ;
sigma_noise = sqrt( noise_power / factor )

disp('Generate random noise n_lmp of std dev sigma')
n_lmp = sigma_noise * randn(P, L^2);
    for en = 1:P
      for el = 0:L-1
        ind = ssht_elm2ind(el, 0);
        n_lmp(en,ind) = sigma_noise .* randn * (en/P) ;
        for m = 1:el
          ind = ssht_elm2ind(el, m);
          n_lmp(en,ind) = ...
            sqrt(sigma_noise^2./2) .* randn * (en/P)  ...
            + sqrt(-1) * sqrt(sigma_noise^2./2) .* randn * (en/P)  ;
        end
      end
    end
    
noise_power = sum(sum(abs(n_lmp).^2));
SNR_before = 10*log10( f_power / noise_power )

disp('Construct the noise in real space...')
initial_noise = flag_synthesis(n_lmp, 'R', R, 'Reality', true);

disp('Add the noise to the initial dataset')
%g_lmp = f_lmp + n_lmp ;
g = f + initial_noise ;

disp('Perform wavelet decomposition...')
[g_wav, g_scal] = flaglet_axisym_analysis(g, 'B_l', B_l, 'B_p', B_p, 'L', L, 'P', P, 'J_min_l', J_min_l, 'J_min_p', J_min_p, 'Reality', true, 'Downsample', Downsample);

disp('Construct noise model...')
noisemodel = cell(J_l+1, J_p+1);
for jl = J_min_l:J_l
    if Downsample == true
        bandlimit_l = min([ ceil(B_l^(jl+1)) L ]);
    else
        bandlimit_l = L;
    end
    for jp = J_min_p:J_p
        if Downsample == true
            bandlimit_p = min([ ceil(B_p^(jp+1)) P ]);
        else
            bandlimit_p = P;
        end
        [jl, jp]
        temp = kappa{jl+1,jp+1};  
        contrib = zeros(1, bandlimit_p);
        nodes = slag_sampling(bandlimit_p, R);
        for r = 0:bandlimit_p-1
            birs = [nodes(r+1), R];
            for n = 0:bandlimit_p-1
                fn = zeros(1,P);
                fn(n+1) = 1.0;
                [K_n_s, ~] = slag_synthesis(fn, 'Nodes', birs) ;
                K_n_s = K_n_s(1);
                for l = 0:bandlimit_l-1
                    contrib(r+1) = ...
                        contrib(r+1) + ((2*l+1)/(4*pi)) * (n/P) * (  K_n_s  * temp(n+1,l+1)).^2 ;
                end
                  
            end
        end
        noisemodel{jl+1,jp+1} = sigma_noise * sqrt(contrib);
    end  
end

g_scal_rec = g_scal ;

disp('Threshold wavelet coefficients...')
g_wav_rec = cell(J_l+1-J_min_l, J_p+1-J_min_p);
nb = 0;
for jl = J_min_l:J_l
    if Downsample == true
        bandlimit_l = min([ ceil(B_l^(jl+1)) L ]);
    else
        bandlimit_l = L;
    end
    for jp = J_min_p:J_p
        if Downsample == true
            bandlimit_p = min([ ceil(B_p^(jp+1)) P ]);
        else
            bandlimit_p = P;
        end
        temp = g_wav{jl+1-J_min_l,jp+1-J_min_p};
        treshold = 3 * noisemodel{jl+1,jp+1};
        [jl, jp]
        for p = 1:2*bandlimit_l-1
            for t = 1:bandlimit_l
                for r = 1:bandlimit_p
                    if abs(temp(r, t, p)) < treshold(r) % Hard thresholding
                        temp(r, t, p) = 0;
                        nb = nb + 1;
                    end
                end
            end
        end
        g_wav_rec{jl+1-J_min_l,jp+1-J_min_p} = temp;
    end
end

disp('Reconstruct the denoised field from the wavelets...')
g_denoised = flaglet_axisym_synthesis(g_wav_rec, g_scal_rec, 'B_l', B_l, 'B_p', B_p, 'L', L, 'P', P, 'J_min_l', J_min_l, 'J_min_p', J_min_p, 'Reality', true, 'Downsample', Downsample);
g_lmp_denoised = flag_analysis(g_denoised, 'Reality', true);

remaining_noise = g_denoised - f;
f_power = sum(sum(abs(f_lmp).^2))

noise_power = sum(sum(abs(n_lmp).^2))
n_lmp_denoised = flag_analysis(remaining_noise, 'Reality', true);
remaining_noise_power = sum(sum(abs(n_lmp_denoised).^2))

SNR_before = 10*log10( f_power / noise_power )
SNR_after = 10*log10( f_power / remaining_noise_power )


Lplot = 192;
Pplot = 192;

% Shell to plot
layer = 173;

disp('Plot...')
figure('Position',[100 100 1000 850]) 

subplot(2,2,1); 
flaglet_plot_f( f, 'L', Lplot, 'P', Pplot, 'layer', layer   )
zoom(zoomfactor)
caxis([minbound maxbound])

subplot(2,2,2);
flaglet_plot_f( initial_noise, 'L', Lplot, 'P', Pplot, 'layer', layer   )
zoom(zoomfactor)
caxis([minbound maxbound])

subplot(2,2,3); 
flaglet_plot_f( g, 'L', Lplot, 'P', Pplot, 'layer', layer   )
zoom(zoomfactor)
caxis([minbound maxbound])

subplot(2,2,4); 
flaglet_plot_f( g_denoised, 'L', Lplot, 'P', Pplot, 'layer', layer   )
zoom(zoomfactor)
caxis([minbound maxbound])

