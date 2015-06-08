function val = var3d( f )

sz = size(f);
d1 = sz(1);
d2 = sz(2);
d3 = sz(3);

vec = [];
for l1 = 0:d1-1
    for l2 = 0:d2-1
        for l3 = 0:d3-1
            vec = [vec f(l1+1,l2+1,l3+1)];
        end
    end
end

val = var(vec);
end