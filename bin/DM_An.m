function f=DM_An(b,h,L,x,y,t,n)
% Calculates An for "nth" acoustic mode using fcs function for Fresnel.
% Takes t as arguments.
    

    epsil = b / L;
    X = (epsil^2)*x;
    Y = epsil*y;
    l = epsil*L;
    ypos = (l + Y) / 2;
    yneg = (l - Y) / 2;
    v = X ./ DM_BigK(h,x,t,n);
    chi = v ./ 2;
    f = abs(((1 - 1i) ./ 2)*(real(fcs(sqrt(2 ./ (pi*chi)) .* ypos)) + real(fcs(sqrt(2 ./ (pi*chi)) .* yneg))) + ((1 + 1i) ./ 2)*(imag(fcs(sqrt(2 ./ (pi*chi)) .* ypos)) + imag(fcs(sqrt(2 ./ (pi*chi)) .* yneg))));
end

    

    