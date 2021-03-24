function [beta0,beta0g] = acoustic_disp_rel_justreal(h,f)
%Ali Abdolali EMC/NCEP/NOAA ali.abdolali@noaa.gov 22, March 2021
%This script caculates the real root of the dispersion relation without 
%gravity term in the governing equation (beta0) and with gravity term in 
%the governing equation beta0g
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see Abdolali and Kirby 2017 for more info
%Abdolali, A., & Kirby, J. T. (2017). Role of compressibility on tsunami 
%propagation. Journal of Geophysical Research: Oceans, 122, 9780? 9794. 
%https://doi.org/10.1002/2017JC013054
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%inputs
%h=water depth (m)
%f frequency (Hz)

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



