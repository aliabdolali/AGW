function f=DM_A0(b,h,L,x,y,t)
% Calculates An for "nth" acoustic mode using fcs function for Fresnel.
% Takes t & n as arguments.
    %global b L
c=1500;
g=10;
     Om_0 = sqrt(2*g/h)*sqrt((sqrt(g*h)*t./x)-1);
    BigK0=Om_0/sqrt(g*h);
    epsil = b / L;
    X = (epsil^2)*x;
    Y = epsil*y;
    l = epsil*L;
    ypos = (l + Y) / 2;
    yneg = (l - Y) / 2;
    v = X ./ BigK0;
    chi = v ./ 2;
    f = abs(((1 - 1i) ./ 2)*(real(fcs(sqrt(2 ./ (pi*chi)) .* ypos)) + real(fcs(sqrt(2 ./ (pi*chi)) .* yneg))) + ((1 + 1i) ./ 2)*(imag(fcs(sqrt(2 ./ (pi*chi)) .* ypos)) + imag(fcs(sqrt(2 ./ (pi*chi)) .* yneg))));
end

    

    