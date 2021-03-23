function [beta0,beta0g,beta_i] = acoustic_disp_rel_gravity(h,f,n_mod)
%Ali Abdolali EMC/NCEP/NOAA ali.abdolali@noaa.gov 22, March 2021
%This script caculates the real root of the dispersion relation without 
%gravity term in the governing equation (beta0); with gravity term in 
%the governing equation (beta0g) and imaginary roots for a given number of
% acoustic modes, n_mod (beta_i) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see Abdolali and Kirby 2017 for more info
%Abdolali, A., & Kirby, J. T. (2017). Role of compressibility on tsunami 
%propagation. Journal of Geophysical Research: Oceans, 122, 9780? 9794. 
%https://doi.org/10.1002/2017JC013054
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%inputs
%h=water depth (m)
%f frequency (Hz)
%n_mod=number of acoustic modes
g=9.81;
eps=0.0000000001;
err=1;
c=1500; %sound speed in water
omega=2*pi*f;
%real root
%classic dispersion
L = g/2/pi*(2*pi/omega)^2*tanh(sinh(omega*(h/g)^0.5));
k0 = 2*pi/(L);
errr=1;
while errr>eps
k0h=k0*h;

ff0=(omega^2)/g/k0-tanh(k0h);
ffp0=-((omega^2/(g*k0^2))+h/(cosh(k0h))^2);
knn0=k0-ff0/ffp0;
errr=abs(k0-knn0)/knn0;
k0=knn0;
end

beta0=k0;

%%%gravity included
%initial guess
k0=omega^2/(g*sqrt(tanh(omega^2*h/g)));
gamma=g/c^2;

errr=1;
while errr>eps
k0h=k0*h;
A=gamma/2/k0;
B=(1-A*tanh(k0h));
C=k0*(1-A^2)*tanh(k0h);
ff0=(omega^2)/g*B-C;
ffp0=tanh(k0h)*((omega^2*A/(g*k0))-(omega^2*h*A*2/(g*sinh(2*k0h)))-1-A^2-(2*k0h*(1-A^2)/(sinh(2*k0h))));

knn0=k0-ff0/ffp0;
errr=abs(k0-knn0)/knn0;
k0=knn0;
end

beta0g=sqrt(k0^2-gamma^2/4);



%%imaginary roots
[beta00,beta_segn] = acoustic_disp_rel(h,f,n_mod);
for i=1:n_mod

ki0=beta_segn(i);


errr=1;


while errr>eps
ki0h=ki0*h;
Ai=gamma/2/ki0;
Bi=(1-Ai*tan(ki0h));
Ci=ki0*(1+Ai^2)*tan(ki0h);
ffi0=g/(omega^2)+Bi/Ci;

ffip0=(tan(ki0h)*(((1-Ai^2)+(2*ki0h*(1+Ai^2)/sin(2*ki0h)))*(-Bi/Ci^2)+((Ai/Bi)-(2*Ai*h/(sin(2*ki0h))))/Ci));

knni0=ki0-ffi0/ffip0;
errr=abs(ki0-knni0)/knni0;
ki0=knni0;
end
beta_i(i)=sqrt(ki0^2+gamma^2/4);
end

