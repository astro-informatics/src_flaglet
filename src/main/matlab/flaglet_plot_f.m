function flaglet_plot_f( f, varargin )

sz = size(f);
N = sz(1);
L = sz(2);

p = inputParser;
p.addParamValue('layer', 0, @isnumeric);
p.addParamValue('L', L, @isnumeric);
p.addParamValue('P', N, @isnumeric);
p.parse(varargin{:});
args = p.Results;

if args.L > L || args.P > N
    nodes = slag_sampling(args.P, 1.0);
    f_ini = f;
    disp('flag analysis')
    f_lmp = flag_analysis(f_ini, 'R', 1.0, 'reality', true);
    f_lmp_zeropadded = zeros(N, args.L^2);
    f_lmp_zeropadded(1:N, 1:L^2) = f_lmp(:,:);
    disp('flag synthesis')
    f = flag_synthesis(f_lmp_zeropadded, 'Nodes', nodes, 'reality', true);
    L = args.L;
    N = args.P ;
end

[rs, thetas, phis] = flag_sampling(L, N, 1.0);
thetas = pi/2 - thetas; 
thetas = [ pi/2 thetas -pi-fliplr(thetas(1:L-1)) pi/2 ];
[rs, thetas] = ndgrid(rs, thetas);
x = rs .* cos(thetas);
y = rs .* sin(thetas);

%ind_neg = find(f < 0);
%f(ind_neg) = f(ind_neg) / abs(min(min(min(f))));
 
%ind_pos = find(f > 0);
%f(ind_pos) = f(ind_pos) / max(max(max(f)));

v_p0 = zeros(N,2*L+1);
v_p0(:,1) = f(:,1,1);
v_p0(:,2:(L+1)) = f(:,:,L);
v_p0(:,(L+2):(2*L)) = fliplr(f(:,1:L-1,1));
v_p0(:,2*L+1) = f(:,1,1);

grid off
hold on
h = surface(x,y,zeros(size(x)),v_p0,'EdgeColor', 'none');
%[C,h] = contourf(x, y, v_p0, nlevels,  'EdgeColor', 'none');
rotate(h,[0 1 0],90)
rotate(h,[1 0 0],90)


v_p90 = zeros(N,2*L+1);
v_p90(:,1) = f(:,1,floor(L/2));
v_p90(:,2:(L+1)) = f(:,:,floor(3*L/2));
v_p90(:,(L+2):(2*L)) = fliplr(f(:,1:L-1,floor(L/2)));
v_p90(:,2*L+1) = f(:,1,floor(L/2));

h = surface(x,y,zeros(size(x)),v_p90,'EdgeColor', 'none');
%[C,h] = contourf(x, y, v_p90, nlevels,  'EdgeColor', 'none');
rotate(h,[0 1 0],90)
rotate(h,[0 0 1],90)
rotate(h,[0 1 0],90)

[rs, thetas, phis] = flag_sampling(L, N, 1.0);
[rs, phis] = ndgrid(rs, [phis phis(1)]);
x = rs .* cos(phis);
y = rs .* sin(phis);
v_t0 = zeros(N,2*L);
for p=1:2*L-1
    v_t0(:,p) = f(:,floor(L/2),p);
end
v_t0(:,2*L) = f(:,floor(L/2),1);

h = surface(x,y,zeros(size(x)),v_t0,'EdgeColor', 'none');
%[C,h] = contourf(x, y, v_t0, nlevels,  'EdgeColor', 'none');
rotate(h,[0 0 1], -90)

view(45,45)
colormap(hot(256))

v = axis;
%axis(0.75*v);
axis off
grid off

rs_line = 1.0;
thetas_line = 0:0.01:2*pi;
[rs_line, thetas_line] = ndgrid(rs_line, thetas_line);
xs_line = rs_line.*cos(thetas_line);
ys_line = rs_line.*sin(thetas_line);
line(xs_line,ys_line,zeros(size(xs_line)),'Color',[0 0 0])
line(zeros(size(xs_line)), xs_line,ys_line,'Color',[0 0 0])
line(xs_line,zeros(size(xs_line)), ys_line,'Color',[0 0 0])


if args.layer ~= 0 && args.layer <= N+1
    layerofinterest = args.layer;
    v_sph = zeros(L, 2*L-1);
    v_sph(:,:) = f(layerofinterest,:,:);

    [rs, thetas, phis] = flag_sampling(L, N, 1.0);
    close = @(x) [x, x(:,1)];
    v_sph = close(v_sph);
    v_norm = (v_sph - min(min(v_sph)))/(max(max(v_sph))- min(min(v_sph)));
    thetas = close(thetas);
    phis = close(phis);

    [thetas, phis] = ndgrid(thetas, phis);
    [x, y, z] = ssht_s2c(thetas, phis, rs(layerofinterest));
    h = surf(x,y,z,v_sph);
    rotate(h,[0 0 1], -90)
    set(h, 'LineStyle', 'none')
    
    rs_line = rs(layerofinterest);
    thetas_line = 0:0.01:2*pi;
    [rs_line, thetas_line] = ndgrid(rs_line, thetas_line);
    xs_line = rs_line.*cos(thetas_line);
    ys_line = rs_line.*sin(thetas_line);
    mycolor = [0 0 0];
    line(xs_line,ys_line,zeros(size(xs_line)),'Color',mycolor)
    line(zeros(size(xs_line)), xs_line,ys_line,'Color',mycolor)
    line(xs_line,zeros(size(xs_line)), ys_line,'Color',mycolor)

    
end

colormap(hsv)

camlight
camlight(-80,10)
material dull

end
