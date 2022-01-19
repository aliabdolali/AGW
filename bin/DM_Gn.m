function f=DM_Gn(b,h,tau,x, t, n)
    T=tau/2;
    f = ( 4 .* sin(DM_BigK(h,x,t,n).*b) .* sin(DM_BigOmega(h,x,t,n).*T) ) ./ (DM_BigK(h,x,t,n) .* DM_BigOmega(h,x,t,n)) ;
  end

    

    