function b3let_plot_axisym_tilling(B_l, B_n, L, N, J_min_l, J_min_n)

% b3let_plot_axisym_tilling - Plot tilling in harmonic space.
% -- Axisymmetric wavelets on the solid sphere.
%
% Default usage :
%
%   b3let_plot_axisym_tilling(B_l, B_n, L, N, J_min_l, J_min_n)
%
% B_l is the wavelet parameter for angular space,
% B_n is the wavelet parameter for radial space,
% L is the angular band-limit,
% N is the radial band-limit,
% J_min_l the first angular wavelet scale to use,
% J_min_n the first radial wavelet scale to use.
%
% B3LET package to perform Wavelet transform on the Solid Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

[kappa kappa0] = b3let_axisym_tilling(B_l, B_n, L, N, J_min_l, J_min_n);
[kappa_l kappa0_l] = s2let_axisym_tilling_mex(B_l, L, J_min_l);
[kappa_n kappa0_n] = s2let_axisym_tilling_mex(B_n, N, J_min_n);

J_l = s2let_jmax(L, B_l);
J_n = s2let_jmax(N, B_n);

figure('Position',[100 100 1000 1000])%figure('Position',[1 1 1000 1000])

colours_n = zeros((J_n+2),3);
hax = axes('Position', [.05, .05, .27, .55]);
colours_n(1,:) = rand(1,3)*0.8+0.2;
plot(0:N-1, kappa0_n, 'k', 'LineWidth', 4, 'Color', colours_n(1,:));
hold on;
for jn = J_min_n:J_n  
  colours_n(jn+2,:) = rand(1,3)*0.8+0.2;
  plot(0:N-1, kappa_n(jn+1,:), 'LineWidth', 4, 'Color', colours_n(jn+2,:));
end
set(gca,'FontSize',20);
axis([0 N-1 0 1.2])
set(gca,'YTick',[0 1]);
set(gca,'LineWidth',4);
%xlabel('x1')
%ylabel('y1')
view(-90,90)
%hold off

colours_l = zeros((J_l+2),3);
hax = axes('Position', [.40, .67, .55, .27]);
colours_l(1,:) = rand(1,3)*0.8+0.2;
plot(0:L-1, kappa0_l, 'k', 'LineWidth', 4, 'Color', colours_l(1,:));
hold on;
for jl = J_min_l:J_l
  colours_l(jl+2,:) = rand(1,3)*0.8+0.2;
  plot(0:L-1, kappa_l(jl+1,:), 'LineWidth', 4, 'Color', colours_l(jl+2,:));
end
set(gca,'FontSize',20);
set(gca,'YTick',[0 1]);
set(gca,'LineWidth',4);
%xlabel('x2')
%ylabel('y2')
axis([0 L-1 0 1.2])
hold off

colours = zeros((J_l+1)*(J_n+1),3);

%figure('Position',[100 100 1000 1000])

hax = axes('Position', [.40, .05, .54, .55]);
colour = 0.5* colours_n(1,:).*colours_l(1,:);
surf(0:L-1, 0:N-1, kappa0, 'FaceColor', colour, 'EdgeColor', 'none');%,'FaceAlpha','flat', ...
        %'AlphaDataMapping','none','AlphaData',(kappa0));
hold on;
for jl = J_min_l:J_l 
    for jn = J_min_n:J_n 
        temp = kappa{jl+1,jn+1};
        colours((jn)*(J_l+1)+jl+1,:) = colours_n(jn+2,:).*colours_l(jl+2,:);%rand(1,3)*0.9;
        surf(0:L-1, 0:N-1, temp, 'FaceColor', colours((jn)*(J_l+1)+jl+1,:), 'EdgeColor', 'none');%,'FaceAlpha','flat', ...
        %'AlphaDataMapping','none','AlphaData',(temp));
    end
end
surf(0:L-1, 0:N-1, zeros(N,L), 'FaceColor', 'black','EdgeColor', 'none')
axis([0 L-1 0 N-1 0 1.2])
%xlabel('x3')
%ylabel('y3')
set(gca,'FontSize',20);
set(gca,'LineWidth',4);
view(0,90)
colormap jet
%alpha(.6)
hold off

hax = axes('Position', [.05, .65, .3, .3]);
surf(0:L-1, 0:N-1, kappa0, 'FaceColor', colour, 'EdgeColor', 'none');%,'FaceAlpha','flat', ...
        %'AlphaDataMapping','none','AlphaData',(kappa0));
hold on;
for jl = J_min_l:J_l 
    for jn = J_min_n:J_n 
        temp = kappa{jl+1,jn+1};
        surf(0:L-1, 0:N-1, temp, 'FaceColor', colours((jn)*(J_l+1)+jl+1,:), 'EdgeColor', 'none');%,'FaceAlpha','flat', ...
        %'AlphaDataMapping','none','AlphaData',(temp));
    end
end
surf(0:L-1, 0:N-1, zeros(N,L), 'FaceColor', 'black','EdgeColor', 'none')
axis([0 L-1 0 N-1 0 1.2])
view(135,45)
set(gca,'FontSize',20);
set(gca,'LineWidth',4);
%xlabel('x4')
%ylabel('y4')
colormap hot(256)
%=alpha(.6)
hold off

end