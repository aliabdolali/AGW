function f=DM_BigOmega(h,x,t,n)
% Calculates stationary phase frequency for "nth" acoustic mode. Takes x, t
% & n as arguments. Outputs "Big Omega (n)".
    %global c h 
    c=1500;
    f = DM_omega(h,n)*( (c.*t) ./ sqrt((c^2).*(t.^2)-(x.^2)) );
end

    

    