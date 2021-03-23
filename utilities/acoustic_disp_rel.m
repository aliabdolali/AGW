function [beta0,beta_segn] = acoustic_disp_rel(h,f,n_mod)
%Ali Abdolali EMC/NCEP/NOAA ali.abdolali@noaa.gov 22, March 2021
%This script caculates the real root of the dispersion relation without 
%gravity term in the governing equation (beta0)and imaginary roots for a 
%given number of acoustic modes, n_mod (beta_i) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see Sammarco et al 2013 for more info
%Sammarco, P., Cecioni, C., Bellotti, G. and Abdolali, A., 2013, Depth-
%integrated equation for large-scale modelling of low-frequency 
%hydroacoustic waves. Journal of Fluid Mechanics, 722, R6 
%doi:10.1017/jfm.2013.153,May 2013.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%inputs
%h=water depth (m)
%f frequency (Hz)
%n_mod=number of acoustic modes
omega=2*pi*f;
gg=9.81;

stampo=0;
if stampo ==1
    figure    
end

% modo zero
L = gg/2/pi*(2*pi/omega)^2*tanh(sinh(omega*(h/gg)^0.5));
beta0 = 2*pi/(L);

% acoustic modes
% n=1:2;
effe = @(x)abs(omega^2*h/9.81/x+tan(x));
for ir=1:n_mod
    effe_aprx = @(x)abs(omega^2*h/9.81/x+1/((n(ir)-1/2)*pi-x));
    dteta=pi/10000;
    kh=(ir-1)*pi+[pi/2:dteta:3/2*pi-dteta];
    left=omega^2.*h/9.81./kh;
    right=-tan(kh);
    [x,fval]=fminbnd(effe,(ir-1)*pi+pi/2,(ir-1)*pi+pi);
    x1 = -omega^2.*h.*pi*(ir-1/2)./(9.81-omega^2.*h); 

        radice(ir)=x;

    if omega^2.*h/9.81/((ir-1)*pi+pi/2) > 5
    radice(ir) = x1;
    end
    if stampo ==1
        semilogy(kh,left)
        hold on
        semilogy(kh,right,'r')
        effep=abs(omega^2*h/9.81./kh+tan(kh));
        effepp=abs(omega^2*h/9.81./kh+(1/(n(ir)-1/2)*pi-x));
        grid on
        semilogy(radice(ir),omega^2*h/9.81/radice(ir),'k.','markersize',15)
    end
end
beta_segn = radice./h;
if stampo == 1
    title(['h=' num2str(h) 'm, f=' num2str(f) 'Hz'])
end