function f = gen_geomodel( L, P )

maxdepth = 3000;

nodes = slag_sampling(P, maxdepth);

for r = 0:P-1
    depth = maxdepth - nodes(r+1);
    cosimod = interpJRmodel(depth,'xpcube');
    [alm_trunc, L_trunc] = cosi2alm( cosimod );
    if r == 0
        f = zeros(P, L, 2*L-1);
    end
    alm = zeros(L^2, 1);
    if L_trunc < L
        alm(1:L_trunc^2, 1) = alm_trunc(:, 1);
    else
        alm(:, 1) = alm_trunc(1:L^2, 1);
    end
    temp = ssht_inverse(alm, L);
    f(r+1,:,:) = real(temp);
end

end