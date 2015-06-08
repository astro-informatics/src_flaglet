% noise_prop_check

niter = 20;

N = 16 ;
L = 16 ;
R = 430.4;

sigma_noise = 0.1 ;

B_l = 2 ;
B_n = 2 ;
J_min_l = 0 ;
J_min_n = 0 ;
J_l = ceil(log(L) ./ log(B_l));
J_n = ceil(log(N) ./ log(B_n));

[kappa kappa0] = b3let_axisym_tilling(B_l, B_n, L, N, J_min_l, J_min_n);

nodes = slag_sampling(N, R);

noisemodel0 = zeros(1,N+1);
for r = 0:N
    birs = [nodes(r+1), R];
    for n = 0:N-1
        fn = zeros(1,N);
        fn(n+1) = 1.0;
        [K_n_s, ~] = slag_synthesis(fn, 'Nodes', birs) ;
        K_n_s = K_n_s(1);
        for l = 0:L-1
            noisemodel0(r+1) = ...
                noisemodel0(r+1) + ((2*l+1)/(4*pi)) * (  K_n_s * kappa0(n+1,l+1)).^2 ;
        end
    end
end
noisemodel0 = sigma_noise * sqrt(noisemodel0);

noisemodel = zeros(J_l+1, J_n+1, N+1);
for jl = 0:J_l
    for jn = 0:J_n
        temp = kappa{jl+1,jn+1};  
        for r = 0:N
            birs = [nodes(r+1), R];
            for n = 0:N-1
                fn = zeros(1,N);
                fn(n+1) = 1.0;
                [K_n_s, ~] = slag_synthesis(fn, 'Nodes', birs) ;
                K_n_s = K_n_s(1);
                for l = 0:L-1
                    noisemodel(jl+1,jn+1,r+1) = ...
                        noisemodel(jl+1,jn+1,r+1) + ((2*l+1)/(4*pi)) * (  K_n_s * temp(n+1,l+1)).^2 ;
                end
                  
            end
        end
    end  
end
noisemodel = sigma_noise * sqrt(noisemodel);

n_scal_std_mean = zeros(1,N+1); 
n_std_mean = zeros(J_l+1, J_n+1, N+1);

for iter = 1:niter
    iter
    
    % Generate random noise n_lmn of std dev sigma
    %n_lmn = sqrt(sigma_noise_ini) * randn(N,L^2);%randreal(N, L);

    n_lmn = sigma_noise * randn(N, L^2);
    for en = 1:N
      for el = 0:L-1
        ind = ssht_elm2ind(el, 0);
        n_lmn(en,ind) = sigma_noise .* randn ;
        for m = 1:el
          ind = ssht_elm2ind(el, m);
          n_lmn(en,ind) = ...
            sqrt(sigma_noise^2./2) .* randn  ...
            + sqrt(-1) * sqrt(sigma_noise^2./2) .* randn ;
        end
      end
    end
    sqrt(var2d(n_lmn));
    %n_lmn = n_lmn * sqrt(sigma_noise_ini/sigma_noise);
    %sigma_noise = var2d(n_lmn);
    noise = flag_synthesis(n_lmn, 'R', R, 'Reality', true);
    [n_wav, n_scal] = b3let_axisym_analysis(noise, 'B_l', B_l, 'B_p', B_n, 'L', L, 'P', N, 'J_min_l', J_min_l, 'J_min_p', J_min_n, 'Reality', true, 'downsample', false);

    % Scaling function
    n_scal_std = zeros(1,N); 
    for r = 0:N
        layer = n_scal(r+1,:,:);
        n_scal_std(r+1) = sqrt(var3d(layer));
    end
    n_scal_std_mean(:) = n_scal_std_mean(:) + n_scal_std(:) ./ niter ;
    
    % Wavelets
    n_std = zeros(N);
    for jl = 0:J_l
        for jn = 0:J_n
            temp = n_wav{jl+1,jn+1};
            for r = 0:N
                layer = temp(r+1,:,:);
                n_std(r+1) = sqrt(var3d(layer));
                n_std_mean(jl+1, jn+1, r+1) = n_std_mean(jl+1, jn+1, r+1) + n_std(r+1) ./ niter ;
            end
           
        end
    end
        
end

noisemodel0 ./ n_scal_std_mean
noisemodel ./ n_std_mean
mean( mean ( mean ( noisemodel ./ n_std_mean ) ) )

