% flaglet_demo2
% Plot one wavelet and its radially translated version

R = 1.0 ;
Rplot1 = 0.49 ;
Rplot2 = 0.49 ;
Rplot3 = 0.49 ;
limit1 = 0.00008;
limit2 = 0.00004 ;
limit3 = 0.000018 ;

s1 = 0.2 ;
s2 = 0.29 ;
s3 = 0.38 ;

L = 92 ;
N = 64 ;
P = 256 ;

B_l = 2 ;
B_n = 2 ;
J_min_l = 0 ;
J_min_n = 0 ;

tau = flag_get_tau(N, R);

[kappa kappa0] = flaglet_axisym_tiling(B_l, B_n, L, N, J_min_l, J_min_n);

% Original
kappa_ln = kappa{5,5};
row = 15;

% First translation
kappa_ln_translated1 = kappa_ln  ;
birs = [ s1, R ];
for n = 1:N
    fn = zeros(1,N) ;
    fn(n) = 1.0 ;
    [K_n_s, ~] = slag_synthesis(fn, 'Nodes', birs) ;
    kappa_ln_translated1(n,:) = kappa_ln(n,:) * K_n_s(1) ;
end

% Second translation
kappa_ln_translated2 = kappa_ln ;
birs = [ s2, R ];
for n = 1:N
    fn = zeros(1,N) ;
    fn(n) = 1.0 ;
    [K_n_s, ~] = slag_synthesis(fn, 'Nodes', birs) ;
    kappa_ln_translated2(n,:) = kappa_ln(n,:) * K_n_s(1) ;
end

% Third translation
kappa_ln_translated3 = kappa_ln ;
birs = [ s3, R ];
for n = 1:N
    fn = zeros(1,N) ;
    fn(n) = 1.0 ;
    [K_n_s, ~] = slag_synthesis(fn, 'Nodes', birs) ;
    kappa_ln_translated3(n,:) = kappa_ln(n,:) * K_n_s(1) ;
end

figure('Position',[100 100 300 1100])

h = (R/P);
nodes = (0:h:R) ;
%nodes = slag_sampling(N, R) ;
%fn_orig = transpose(kappa_ln(:,row)) ;
%[f ~] = slag_synthesis(fn_orig, 'Nodes', nodes);
%subplot(3,2,2)
%plot(nodes, f, 'black')
%xlabel('Radius')
%ylabel('Amplitude')
%title('Radial profile')
%axis([0.01 Rplot1 -limit1 limit1])
%[ min(f) max(f) mean(f) ]

fn_transl1 = transpose(kappa_ln_translated1(:,row)) ;
[f_transl1 ~] = slag_synthesis(fn_transl1, 'Nodes', nodes);
subplot(3,1,1)
plot(nodes, f_transl1, 'black')
xlabel('x01')
ylabel('y01')
axis([0.01 Rplot1 -limit1 limit1])
%[ min(f_transl1) max(f_transl1) mean(f_transl1) ]


fn_transl2 = transpose(kappa_ln_translated2(:,row)) ;
[f_transl2 ~] = slag_synthesis(fn_transl2, 'Nodes', nodes);
subplot(3,1,2)
plot(nodes, f_transl2, 'black')
xlabel('x02')
ylabel('y02')
axis([0.01 Rplot2 -limit2 limit2])
%[ min(f_transl2) max(f_transl2) mean(f_transl2) ]

fn_transl3 = transpose(kappa_ln_translated3(:,row)) ;
[f_transl3 ~] = slag_synthesis(fn_transl3, 'Nodes', nodes);
subplot(3,1,3)
plot(nodes, f_transl3, 'black')
xlabel('x03')
ylabel('y03')
axis([0.01 Rplot3 -limit3 limit3])
%[ min(f_transl2) max(f_transl2) mean(f_transl2) ]

colormap(hot(256))

%subplot(3,2,1)
%flaglet_plot_axisym_wavelet_kernel(kappa_ln, R, Rplot1, P)
%title('Original wavelet')

figure('Position',[100 100 300 1100])

subplot(3,1,1)
flaglet_plot_axisym_wavelet_kernel(kappa_ln_translated1, R, Rplot1, P)
%title(['Wavelet translated by ',num2str(s1)])

subplot(3,1,2)
flaglet_plot_axisym_wavelet_kernel(kappa_ln_translated2, R, Rplot2, P)
%title(['Wavelet translated by ',num2str(s2)])

subplot(3,1,3)
flaglet_plot_axisym_wavelet_kernel(kappa_ln_translated3, R, Rplot3, P)
%title(['Wavelet translated by ',num2str(s3)])

%colormap(CubeHelix(256,0.5,-1.5,1.2,1.0))
colormap(hot(256))