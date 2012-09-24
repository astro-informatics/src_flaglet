function flaglet_plot_axisym_wavelet_kernel(kappa_ln, R, Rplot, P)

sz = size(kappa_ln);
N = sz(1);
L = sz(2);

flmn = zeros(N,L^2);
for l = 0:L-1
    for n = 0:N-1
        flmn(n + 1, l^2 + l + 1) = kappa_ln(n+1,l+1);
    end
end

h = (R/(P));
rs = 0:h:R ;
[~, thetas, ~] = flag_sampling(L, N, 1.0);

ind = find(rs <= Rplot);
nind = length(ind);

thetas = pi/2 - thetas; 
thetas = [ pi/2 thetas -pi-fliplr(thetas(1:L-1)) pi/2 ];

f = (flag_synthesis(flmn, 'Nodes', rs, 'reality', true));
f(1:3,:,:) = 0.0;

[rs, thetas] = ndgrid(rs(ind),thetas);
x = rs .* cos(thetas);
y = rs .* sin(thetas);

v = zeros(nind,2*L+1);
v(:,1) = f(ind,1,1);
v(:,2:(L+1)) = f(ind,:,L);
v(:,(L+2):(2*L)) = fliplr(f(ind,1:L-1,1));
v(:,2*L+1) = f(ind,1,1);

%size(x)
%size(y)
%size(v)

ind_neg = find(v < 0);
v(ind_neg) = v(ind_neg) / abs(min(min(min(v))));
 
ind_pos = find(v > 0);
v(ind_pos) = v(ind_pos) / max(max(max(v)));

%[min(min(min(v))), max(max(max(v)))] 
%limit = max( [ -min(min(min(v))) max(max(max(v))) ]) ;
%levels = -limit:(2*limit/nlevels):limit ;
%levels = -1:(2/nlevels):1 ;

p = polar([0 2*pi], [0 Rplot]);
grid off
hold on
surface(x, y, zeros(size(x)), v, 'EdgeColor', 'none');
%colorbar
ph=findall(gca,'type','text'); 
ps=get(ph,'string'); 
%disp([num2cell(1:numel(ps)).',ps]); 
ps(1:numel(ps))={''}; 
%ps(16:19)={''};
ps(4:11)={ 
          '90'
          '-90'
          ''
          ''
          ''
          ''
          '180'
          '0'
     }; 

set(ph,{'string'},ps); 
delete(findall(ancestor(p,'figure'),'HandleVisibility','off','type','line'));
delete(p)
%axis([-Rplot Rplot -Rplot+0.2 Rplot+0.2])
end