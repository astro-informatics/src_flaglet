% flaglet_perf_tests

ell = [ 4 8 16 32 64 128 256 512 ];
accuracy = [ ...
    3.49e-15 ... % 4
    2.64e-14 ... % 8
    9.02e-14 ... % 16
    9.08e-13 ... % 32
    2.93e-12 ... % 64
    2.56e-11 ... % 128
    1.1e-10 ... % 256
    4.56e-10 ... % 512
    ];
speed = [ ...
    0.0036 ...  % 4
    0.022 ... % 8
    0.3 ... % 16
    3.3 ... % 32
    50.0 ... % 64
    880 ... % 128
    10000 ... % 256
    110000 ... % 512
    ];
accuracy_multires = [ ...
    3.31e-15 ... % 4
    2.27e-14 ... % 8
    8.42e-14 ... % 16
    8.15e-13 ... % 32
    2.86e-12 ... % 64
    2.35e-11 ... % 128
    1.1e-10 ... % 256
    4.56e-10 ... % 512
    ];
speed_multires = [ ...
    0.00030 ...  % 4
    0.005 ... % 8
    0.08 ... % 16
    0.65 ... % 32
    8.0 ... % 64
    120 ... % 128
    2300 ... % 256
    23000 ... % 512
    ];

nb_samples = 2*ell.^3;
nb_samples_pow = str2mat('t04', 't07', 't10', 't13', 't16', 't19', 't22', 't25', 't28');

figure('Position',[100 100 700 700])

subplot(2,1,1)
loglog( nb_samples, (1e-14)*ell.^2, 'red','LineWidth', 2 )
hold on
loglog( nb_samples, accuracy, '--oblack', 'LineWidth', 2, 'MarkerSize',10)
loglog( nb_samples, accuracy_multires, '-oblack', 'LineWidth', 2, 'MarkerSize',10)
grid on
xlabel('x1','FontSize',20)
ylabel('y1','FontSize',20)
set(gca,'XTick',nb_samples,'XTickLabel',nb_samples_pow,'FontSize',20)
set(gca,'YTick',[1e-16 1e-14 1e-12  1e-10 1e-8] ,'FontSize',20)
axis([80 1.5*nb_samples(8) 1e-16 1e-8])

%title('Accuracy of overall transform','FontSize',20)

subplot(2,1,2)
loglog( nb_samples, (5e-5)*ell.^4, 'red','LineWidth', 2 )
hold on
loglog( nb_samples, speed, '--oblack', 'LineWidth', 2, 'MarkerSize',10)
loglog( nb_samples, speed_multires, '-oblack', 'LineWidth', 2, 'MarkerSize',10)
grid on
xlabel('x2','FontSize',20)
ylabel('y2','FontSize',20)
axis([80 1.5*nb_samples(8) 1e-5 2e7])
set(gca,'YTick',[1e-6 1e-4 1e-2  1e0 1e2 1e4 1e6] ,'FontSize',20)
set(gca,'XTick',nb_samples,'XTickLabel',nb_samples_pow,'FontSize',20)
%title('Speed of overall transform','FontSize',20)

