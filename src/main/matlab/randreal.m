function flmn = randreal(N, L)

flmn = zeros(N, L^2);
flmn = randn(size(flmn)) + sqrt(-1)*randn(size(flmn));
flmn = 2.*(flmn - (1+sqrt(-1))./2);
% Impose reality on flms.
for en = 1:N
   for el = 0:L-1
      ind = el*el + el + 1;
      flmn(en,ind) = real(flmn(en,ind));
      for m = 1:el
         ind_pm = el*el + el + m + 1;
         ind_nm = el*el + el - m + 1;
         flmn(en,ind_nm) = (-1)^m * conj(flmn(en,ind_pm));
      end  
   end
end

end