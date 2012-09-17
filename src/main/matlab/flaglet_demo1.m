% flaglet_demo1
% Plot wavelets functions

s = 0.3;

R = 1.0 ;
Rplot = 0.39 ;

L = 92 ;
N = 92 ;
P = 128 ;

scales_l = [4 5];
scales_n = [5 6];

J_l = numel(scales_l);
J_n = numel(scales_n);


B_l = 2 ;
B_n = 2 ;
J_min_l = 0 ;
J_min_n = 0 ;


% Translation
K_n_s = zeros(N,L)  ;
birs = [ s, R ];
for n = 1:N
    fn = zeros(1,N) ;
    fn(n) = 1.0 ;
    [K_n_s_temp, ~] = slag_synthesis(fn, 'Nodes', birs) ;
    K_n_s(n,:) = K_n_s_temp(1) ;
end

[kappa kappa0] = flaglet_axisym_tiling(B_l, B_n, L, N, J_min_l, J_min_n);

figure('Position',[100 100 750 800])
for jl=scales_l
    for jn=scales_n
        ind = (find(jl == scales_l)-1)*(J_n) + find(jn == scales_n);
        subplot(J_l, J_n, ind)
        temp = kappa{jl,jn} ;
        temp = temp .* K_n_s ; 
        flaglet_plot_axisym_wavelet_kernel(temp, R, Rplot, P)
        %title(['jl=',int2str(jl),' jn=',int2str(jn)])
        caxis([-1 1])
        %xlabel(['x',int2str(ind)])
    end
end
colormap(hot(256))


