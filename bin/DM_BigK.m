function f=DM_BigK(h,x,t,n)
% Calculates Kn for the "nth" acoustic mode. Takes x, t & n as arguments.
  %  global c h
  c=1500;
    f = (DM_omega(h,n).*x) ./ ( (c.^2).*t .* sqrt(1-(x./(c.*t)).^2)  );
end

    

    