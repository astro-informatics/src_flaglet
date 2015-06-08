% flaglet_demo5
% Plot wavelets functions

s = 0.15;
Rplot = 0.20 ;
zoomfactor = 2.1;
R = 1.0 ;

L = 64 ;
P = 128 ;
N = 5;
spin = 2;
clayer = 85
layer = clayer;
layers = 80:90;

oversampfac = 2;

B_l = 3 ;
B_n = 3 ;

if true
    J_min_l = 1 ;
    J_min_n = 2 ;
    J_l = ceil(log(L)/log(B_l));
    J_n = ceil(log(P)/log(B_n));
    scales_l = [2 3 4]%J_min_l+1:J_l+1;
    scales_n = [3 4]%J_min_n+1:J_n+1
else
    scales_l = [2 3];
    scales_n = [3 4];
    J_l = numel(scales_l);
    J_n = numel(scales_n);
end

pltroot = '../../../figs'
configstr = ['Spin',int2str(spin),'_N',int2str(N),'_L',int2str(L),'_Bl',int2str(B_l),'_P',int2str(P),'_Bp',int2str(B_n)]



% Translation
K_n_s = zeros(1,P)  ;
K_ln_s = zeros(P,L*L)  ;
birs = [ s, R ];
for n = 1:P
    fn = zeros(1,P) ;
    fn(n) = 1.0 ;
    [K_n_s_temp, ~] = slag_synthesis(fn, 'Nodes', birs) ;
    K_n_s(1, n) = K_n_s_temp(1) ;
    K_ln_s(n,:) = K_n_s_temp(1) ;
end

[kappa kappa0] = flaglet_tiling(B_l, B_n, L, P, N, spin, J_min_l, J_min_n);

h = (1.0/(oversampfac*P));
nodes = (0:h:1.0) ;
taufac = flag_get_tau(P, R) / flag_get_tau(oversampfac*P, R)
nodes = slag_sampling(oversampfac*P, Rplot);

nx = 4;
ny = 2;


fig = figure('Renderer','zbuffer', 'Position',[0 0 1000 550],'Color',[1 1 1]);
%set(gca,'NextPlot','replaceChildren');
%file = [pltroot,'/flaglet_demo5_', configstr, '_wav_hot_movie.png']
%writerObj = VideoWriter( file );
%writerObj.FileFormat
%writerObj.FrameRate = 3;
%writerObj.Quality = 100;
%open(writerObj);


reality = spin == 0


%nblayers = length(layers)
%for ilayer = 1:nblayers
    %layer = layers(ilayer)
    indfig = 1
    for jl=scales_l
        for jn=scales_n
            if indfig <= nx * ny    
                subplot(ny, nx, indfig)
                f = zeros(oversampfac*P,L,(2*L-1));
                f_lmp = kappa{jl,jn} .* K_ln_s(:,:);
                f_lmp_over = zeros(oversampfac*P,L*L);
                f_lmp_over(1:P,:) = f_lmp(:,:);
                [jl, jn]
                bounds_kappa = [min(min(min(kappa{jl,jn}))), max(max(max(kappa{jl,jn})))]
                bounds_f_lmp = [min(min(min(f_lmp))), max(max(max(f_lmp)))]
                if spin == 0
                    f = flag_synthesis(f_lmp_over, 'reality', true, 'Nodes', nodes);%
                    bounds_f = [min(min(min(f))), max(max(max(f)))]
                    flaglet_plot_f(f, 'ShowSlices', true, 'Rplot', Rplot / taufac, 'layer', layer) 
                    %title(['Wavelet : jl=',int2str(jl),' jn=',int2str(jn)])
                    %if ilayer == 1
                        colormap(hot)   
                        caxis([min(min(f(clayer,:,:))), max(max(f(clayer,:,:)))]) 
                        zoom(zoomfactor)
                        view(45,45)
                    %end
                    indfig = indfig + 1;
                else
                    freal = zeros(oversampfac*P,L,(2*L-1));
                    fimag = zeros(oversampfac*P,L,(2*L-1));
                    for p = 1:P
                        temp = ssht_inverse(f_lmp_over(p,:), L, 'Spin', spin);
                        for p1 = 1:L
                            for p2 = 1:2*L-1
                                freal(p,p1,p2) = real(temp(p1,p2));
                                fimag(p,p1,p2) = imag(temp(p1,p2));
                            end
                        end
                    end
                    for p1 = 1:L
                        for p2 = 1:2*L-1
                            [tempreal, nodes] = slag_synthesis(freal(:,p1,p2), 'P', oversampfac*P, 'Nodes', nodes);
                            [tempimag, nodes] = slag_synthesis(fimag(:,p1,p2), 'P', oversampfac*P, 'Nodes', nodes);
                            for p = 1:oversampfac*P
                                freal(p,p1,p2) = tempreal(p);
                                fimag(p,p1,p2) = tempimag(p);
                            end
                        end
                    end
                    flaglet_plot_f(freal, 'ShowSlices', true, 'Rplot', Rplot / taufac, 'layer', layer) 
                    caxis([min(min(freal(layer,:,:))), max(max(freal(layer,:,:)))]) 
                    zoom(zoomfactor)
                    view(45,45)
                    indfig = indfig + 1
                    subplot(ny, nx, indfig)
                    flaglet_plot_f(fimag, 'ShowSlices', true, 'Rplot', Rplot / taufac, 'layer', layer) 
                    caxis([min(min(fimag(layer,:,:))), max(max(fimag(layer,:,:)))]) 
                    zoom(zoomfactor)
                    view(45,45)
                    indfig = indfig + 1
                end 
            end 
        end
    end
    %F(ilayer) = getframe(fig);
    %frame = getframe(fig);
    %writeVideo(writerObj,frame); 
%end
%close(writerObj);


colormap(hot)
fname = [pltroot,'/flaglet_demo5_', configstr, '_wav_hot.png']
print('-r200', '-dpng', fname)

colormap(jet)
fname = [pltroot,'/flaglet_demo5_', configstr, '_wav_jet.png']
print('-r200', '-dpng', fname)


stop

fields = kappa;
for jl=scales_l
    for jn=scales_n
        temp = kappa{jl,jn} ;
        temp = temp .* K_n_s ;
        flmn = zeros(P,L^2);
        for l = 0:L-1
            for n = 0:P-1
                flmn(n + 1, l^2 + l + 1) = temp(n+1,l+1);
            end
        end
        v = flag_synthesis(flmn,'reality', true);%
        ind_neg = find(v < 0);
        v(ind_neg) = v(ind_neg) / abs(min(min(min(v))));
        ind_pos = find(v > 0);
        v(ind_pos) = v(ind_pos) / max(max(max(v)));
        fields{jl,jn} = v;
    end
end

layer = round(2.99*P/5)
fig = figure('Position',[0 0 1600 900],'Color',[1 1 1]);
for jl=scales_l
    for jn=scales_n
        ind = (find(jn == scales_n)-1)*(J_l) + find(jl == scales_l);
        subplot(J_n, J_l, ind)
        flaglet_plot_f(fields{jl,jn}, 'layer', layer, 'Rplot', Rplot, 'ShowSlices', true);
        colormap(hot)
        zoom(2.1)
        view(45,45)
light;
lighting phong;
    end
end

print -r200 -dtiff flaglet_demo5_highdefplot1.tif
colormap(jet)
print -r200 -dtiff flaglet_demo5_highdefplot2.tif

stop



fig = figure('Renderer','zbuffer', 'Position',[0 0 1000 550],'Color',[1 1 1]);
file = ['waveletsL',int2str(L),'P',int2str(P)]
set(gca,'NextPlot','replaceChildren');
writerObj = VideoWriter( file );
writerObj.FileFormat
writerObj.FrameRate = 3;
writerObj.Quality = 100;
open(writerObj);

nblayers = length(layers)
for ilayer = 1:nblayers
    layer = layers(ilayer)
    for jl=scales_l
        for jn=scales_n
            ind = (find(jn == scales_n)-1)*(J_l) + find(jl == scales_l);
            subplot(J_n, J_l, ind)
            flaglet_plot_f(fields{jl,jn}, 'layer', layer, 'Rplot', Rplot, 'ShowSlices', false);
            if ilayer == 1
                colormap(hot)
                zoom(zoomfactor)
                caxis([-1, 1])
            end
        end
    end
    F(ilayer) = getframe(fig);
    frame = getframe(fig);
    writeVideo(writerObj,frame); 
end
close(writerObj);


fig = figure('Renderer','zbuffer', 'Position',[0 0 1000 550],'Color',[1 1 1]);
file = ['waveletsL',int2str(L),'P',int2str(P), 'ShowSlices']
set(gca,'NextPlot','replaceChildren');
writerObj = VideoWriter( file );
writerObj.FileFormat
writerObj.FrameRate = 3;
writerObj.Quality = 100;
open(writerObj);
nblayers = length(layers)
for ilayer = 1:nblayers
    layer = layers(ilayer)
    for jl=scales_l
        for jn=scales_n
            ind = (find(jn == scales_n)-1)*(J_l) + find(jl == scales_l);
            subplot(J_n, J_l, ind)
            flaglet_plot_f(fields{jl,jn}, 'layer', layer, 'Rplot', Rplot, 'ShowSlices', true);
            if ilayer == 1
                colormap(hot)
                zoom(zoomfactor)
                caxis([-1, 1])
            end
        end
    end
    F(ilayer) = getframe(fig);
    frame = getframe(fig);
    writeVideo(writerObj,frame); 
end
close(writerObj);


% 
% for jl=scales_l
%     for jn=scales_n
%         ind = (find(jl == scales_l)-1)*(J_n) + find(jn == scales_n);
%         %subplot(J_l, J_n, ind)
%         temp = kappa{jl,jn} ;
%         temp = temp .* K_n_s ; 
%         flmn = zeros(P,L^2);
%         for l = 0:L-1
%             for n = 0:P-1
%                 flmn(n + 1, l^2 + l + 1) = temp(n+1,l+1);
%             end
%         end
%         v = flag_synthesis(flmn,'reality', true);% 
%         ind_neg = find(v < 0);
%         v(ind_neg) = v(ind_neg) / abs(min(min(min(v))));
%         ind_pos = find(v > 0);
%         v(ind_pos) = v(ind_pos) / max(max(max(v)));
%         file = ['waveletL',int2str(L),'P'+int2str(P),'wavjl',int2str(jl),'wavjp',int2str(jn)]
%         flagmovie( v, file, L, P, layers, zoomfactor, caxisvals );
%         %flaglet_plot_f(v);
%         %caxis([-1 1]);
%         %zoom(zoomfactor)
%     end
% end
% 
% 
