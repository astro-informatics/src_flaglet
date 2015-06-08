function plot_delta

P = 256
R = 1.0
R_plot = 1.0
rs = [0.6, 0.8, 0.9]
npoints = 5000;
finegrid = 0:R/npoints:R;
colours = rand(3,3)*0.9;
tau = flag_get_tau(P, R);

figure('Position',[100 100 800 200]) 
% subplot(2,1,1)
hold on   
for r = rs
    Kpr = slag_basis(P, r, tau);
    %[delta, nodes] = slag_synthesis(Kpr);
    [delta, ~] = slag_synthesis(Kpr, 'nodes', finegrid);
    plot(finegrid, delta, 'Color', colours(r == rs,:), 'LineWidth', 2)
    %axis([0.1 R_plot -30*10^-7 100*10^-7])
    if r == rs(1) 
        lim = 1.1*max(delta(finegrid > 0.1))
        axis([0.1 R_plot -lim lim]) 
    end
    set(gca, 'box','on')
end
% subplot(2,1,2)
% hold on   
% for r = rs
%     Kpr = slag_basis(P, r, tau);
%     %[delta, nodes] = slag_synthesis(Kpr);
%     [delta, ~] = slag_synthesis(Kpr, 'nodes', finegrid);
%     plot(finegrid, (finegrid).*delta, 'Color', colours(r == rs,:), 'LineWidth', 2)
%     axis([0.1 R_plot -1*10^-5 3.5*10^-5])  
%     set(gca, 'box','on')
% end

end