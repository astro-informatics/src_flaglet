function val = var2d( f )

sz = size(f);
d1 = sz(1);
d2 = sz(2);

vec = [];
for l1 = 0:d1-1
    for l2 = 0:d2-1
        vec = [vec f(l1+1,l2+1)];
    end
end

val = var(vec);
end