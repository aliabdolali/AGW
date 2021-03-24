function [beta_i,q_i,s_i,k_i,c_i,val] = acoustic_disp_rel_elasticbed_imaginaryPP(h,f,lambda,mu,rhol,rhos,Cl,n_mod,plott)
%Ali Abdolali EMC/NCEP/NOAA ali.abdolali@noaa.gov 22, March 2021
%This script caculates the imaginary root of the dispersion relation without 
%gravity term in the governing equation for a given number of
% acoustic modes, n_mod (beta_i) for the case of elastic halfspace sea bottom 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see Abdolali et al 2019 for more info
%Abdolali, A., Kadri, U. & J.T. Kirby,  2019, Effect of Water
%Compressibility, Sea-floor Elasticity, and Field Gravitational Potential 
%on Tsunami Phase Speed, Scientific Reports, Nature, 
%doi:10.1038/s41598-019-52475-0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%inputs
%h=water depth (m)
%f frequency (Hz)
%n_mod=number of acoustic modes
%inputs example
%rhos=2750;%kg/m^3
%rhol=1000;%kg/m^3
%Cl=1500;%m/s %sound speed in water
%Cs=3550;%m/s
%mu=Cs^2*rhos;
%lambda=3.983375000000000e+10;
%Cp=sqrt((lambda+2*mu)/rhos);
%n_mod=4;
%h=4000;
%f=.1;
%plott=1 to see dispersion relation graphyically
%plott=0 to skip graphyical dispersion relation

g=9.81;
omega=2*pi*f;
Cs=sqrt(mu/rhos);%m/s
Cp=sqrt((lambda+2*mu)/rhos);
eps=0.0000000001;
err=1;

omega=2*pi*f;


%%%%%%%%%%%%%%imaginary roots%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%imaginary roots%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%imaginary roots%%%%%%%%%%%%%%%%%%%%%%%%
%initial guess
RR0=pi/2/h;
DRR=pi/h;
R0=RR0:DRR:(2*n_mod-1)*pi/2/h;
R0(2:end+1)=R0;
R0(1)=pi/64/h;
ff = @(x)abs(-tanh((1i*x)*h)+((omega^2/(1i*x))*((sqrt((1i*x)^2+omega^2/Cl^2-omega^2/Cp^2)*rhol*((((1i*x)^2+omega^2/Cl^2)-((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2))/(((1i*x)^2+omega^2/Cl^2)+((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2))))+((1/g)*(((4*mu*((1i*x)^2+omega^2/Cl^2)*sqrt((1i*x)^2+omega^2/Cl^2-omega^2/Cp^2)*sqrt((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2))/(((1i*x)^2+omega^2/Cl^2)+((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2)))-(((lambda+2*mu)*((1i*x)^2+omega^2/Cl^2-omega^2/Cp^2))-(lambda*((1i*x)^2+omega^2/Cl^2)))))))/(((rhol*omega^4*sqrt((1i*x)^2+omega^2/Cl^2-omega^2/Cp^2))/(g*(1i*x)^2))*((((1i*x)^2+omega^2/Cl^2)-((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2))/(((1i*x)^2+omega^2/Cl^2)+((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2)))+(((4*mu*((1i*x)^2+omega^2/Cl^2)*sqrt((1i*x)^2+omega^2/Cl^2-omega^2/Cp^2)*sqrt((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2))/(((1i*x)^2+omega^2/Cl^2)+((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2)))-(((lambda+2*mu)*((1i*x)^2+omega^2/Cl^2-omega^2/Cp^2))-(lambda*((1i*x)^2+omega^2/Cl^2))))));
fff = @(x)(-tanh((1i*x)*h)+((omega^2/(1i*x))*((sqrt((1i*x)^2+omega^2/Cl^2-omega^2/Cp^2)*rhol*((((1i*x)^2+omega^2/Cl^2)-((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2))/(((1i*x)^2+omega^2/Cl^2)+((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2))))+((1/g)*(((4*mu*((1i*x)^2+omega^2/Cl^2)*sqrt((1i*x)^2+omega^2/Cl^2-omega^2/Cp^2)*sqrt((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2))/(((1i*x)^2+omega^2/Cl^2)+((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2)))-(((lambda+2*mu)*((1i*x)^2+omega^2/Cl^2-omega^2/Cp^2))-(lambda*((1i*x)^2+omega^2/Cl^2)))))))/(((rhol*omega^4*sqrt((1i*x)^2+omega^2/Cl^2-omega^2/Cp^2))/(g*(1i*x)^2))*((((1i*x)^2+omega^2/Cl^2)-((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2))/(((1i*x)^2+omega^2/Cl^2)+((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2)))+(((4*mu*((1i*x)^2+omega^2/Cl^2)*sqrt((1i*x)^2+omega^2/Cl^2-omega^2/Cp^2)*sqrt((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2))/(((1i*x)^2+omega^2/Cl^2)+((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2)))-(((lambda+2*mu)*((1i*x)^2+omega^2/Cl^2-omega^2/Cp^2))-(lambda*((1i*x)^2+omega^2/Cl^2))))));
fffRight = @(x)(((omega^2/(1i*x))*((sqrt((1i*x)^2+omega^2/Cl^2-omega^2/Cp^2)*rhol*((((1i*x)^2+omega^2/Cl^2)-((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2))/(((1i*x)^2+omega^2/Cl^2)+((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2))))+((1/g)*(((4*mu*((1i*x)^2+omega^2/Cl^2)*sqrt((1i*x)^2+omega^2/Cl^2-omega^2/Cp^2)*sqrt((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2))/(((1i*x)^2+omega^2/Cl^2)+((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2)))-(((lambda+2*mu)*((1i*x)^2+omega^2/Cl^2-omega^2/Cp^2))-(lambda*((1i*x)^2+omega^2/Cl^2)))))))/(((rhol*omega^4*sqrt((1i*x)^2+omega^2/Cl^2-omega^2/Cp^2))/(g*(1i*x)^2))*((((1i*x)^2+omega^2/Cl^2)-((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2))/(((1i*x)^2+omega^2/Cl^2)+((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2)))+(((4*mu*((1i*x)^2+omega^2/Cl^2)*sqrt((1i*x)^2+omega^2/Cl^2-omega^2/Cp^2)*sqrt((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2))/(((1i*x)^2+omega^2/Cl^2)+((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2)))-(((lambda+2*mu)*((1i*x)^2+omega^2/Cl^2-omega^2/Cp^2))-(lambda*((1i*x)^2+omega^2/Cl^2))))));
fffLeft = @(x)(tanh((1i*x)*h));

%%
for j=1:length(R0)-1
xx=R0(j):(R0(j+1)-R0(j))/20000:R0(j+1);
for l=1:length(xx)
funn(l)=ff(xx(l));
end

[im,jm]=find(funn==min(funn));

R(j)=xx(jm);
funnn(j)=ff(R(j));
end
k0=sqrt(-R.^2+omega^2/Cl^2);
q0=sqrt(k0.^2-omega^2/Cp^2);
s0=sqrt(k0.^2-omega^2/Cs^2);
n_effective=1;
for t=2:n_mod
    if imag(s0(t))==0;
     n_effective=n_effective+1;   
    end
end
funnn;
beta_i=nan*ones(1,n_mod);
val=nan*ones(1,n_mod);
k_i=nan*ones(1,n_mod);
q_i=nan*ones(1,n_mod);
s_i=nan*ones(1,n_mod);
c_i=nan*ones(1,n_mod);

%%

 for i=1:n_effective
    x=R(i);

errr=1;
while errr>eps
%k=sqrt(-r^2+omega^2/Cl^2);
%q=sqrt(-r^2+omega^2/Cl^2-omega^2/Cp^2);
%s=sqrt(-r^2+omega^2/Cl^2-omega^2/Cs^2);


ff0=(-tanh((1i*x)*h)+((omega^2/(1i*x))*((sqrt((1i*x)^2+omega^2/Cl^2-omega^2/Cp^2)*rhol*((((1i*x)^2+omega^2/Cl^2)-((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2))/(((1i*x)^2+omega^2/Cl^2)+((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2))))+((1/g)*(((4*mu*((1i*x)^2+omega^2/Cl^2)*sqrt((1i*x)^2+omega^2/Cl^2-omega^2/Cp^2)*sqrt((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2))/(((1i*x)^2+omega^2/Cl^2)+((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2)))-(((lambda+2*mu)*((1i*x)^2+omega^2/Cl^2-omega^2/Cp^2))-(lambda*((1i*x)^2+omega^2/Cl^2)))))))/(((rhol*omega^4*sqrt((1i*x)^2+omega^2/Cl^2-omega^2/Cp^2))/(g*(1i*x)^2))*((((1i*x)^2+omega^2/Cl^2)-((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2))/(((1i*x)^2+omega^2/Cl^2)+((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2)))+(((4*mu*((1i*x)^2+omega^2/Cl^2)*sqrt((1i*x)^2+omega^2/Cl^2-omega^2/Cp^2)*sqrt((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2))/(((1i*x)^2+omega^2/Cl^2)+((1i*x)^2+omega^2/Cl^2-omega^2/Cs^2)))-(((lambda+2*mu)*((1i*x)^2+omega^2/Cl^2-omega^2/Cp^2))-(lambda*((1i*x)^2+omega^2/Cl^2))))));


ffp0 =-((1i*omega^2*(-((omega^2*rhol*x)/(Cs^2*((2*omega^2)/Cl^2 - omega^2/Cs^2 - 2*x^2)*...
         sqrt(omega^2/Cl^2 - omega^2/Cp^2 - x^2))) + ...
      (4*omega^2*rhol*x*sqrt(omega^2/Cl^2 - omega^2/Cp^2 - x^2))/...
       (Cs^2*((2*omega^2)/Cl^2 - omega^2/Cs^2 - 2*x^2)^2) + ...
      (1/g)*(-2*lambda*x + 2*(lambda + 2*mu)*x - (4*mu*x*(omega^2/Cl^2 - x^2)*...
          sqrt(omega^2/Cl^2 - omega^2/Cp^2 - x^2))/(((2*omega^2)/Cl^2 - omega^2/Cs^2 - ...
           2*x^2)*sqrt(omega^2/Cl^2 - omega^2/Cs^2 - x^2)) - ...
        (4*mu*x*(omega^2/Cl^2 - x^2)*sqrt(omega^2/Cl^2 - omega^2/Cs^2 - x^2))/...
         (((2*omega^2)/Cl^2 - omega^2/Cs^2 - 2*x^2)*sqrt(omega^2/Cl^2 - omega^2/Cp^2 - ...
            x^2)) - (8*mu*x*sqrt(omega^2/Cl^2 - omega^2/Cp^2 - x^2)*...
          sqrt(omega^2/Cl^2 - omega^2/Cs^2 - x^2))/((2*omega^2)/Cl^2 - omega^2/Cs^2 - ...
          2*x^2) + (16*mu*x*(omega^2/Cl^2 - x^2)*sqrt(omega^2/Cl^2 - omega^2/Cp^2 - x^2)*...
          sqrt(omega^2/Cl^2 - omega^2/Cs^2 - x^2))/((2*omega^2)/Cl^2 - omega^2/Cs^2 - ...
           2*x^2)^2)))/(x*(lambda*(omega^2/Cl^2 - x^2) - ...
      (omega^6*rhol*sqrt(omega^2/Cl^2 - omega^2/Cp^2 - x^2))/...
       (Cs^2*g*x^2*((2*omega^2)/Cl^2 - omega^2/Cs^2 - 2*x^2)) - ...
      (lambda + 2*mu)*(omega^2/Cl^2 - omega^2/Cp^2 - x^2) + ...
      (4*mu*(omega^2/Cl^2 - x^2)*sqrt(omega^2/Cl^2 - omega^2/Cp^2 - x^2)*...
        sqrt(omega^2/Cl^2 - omega^2/Cs^2 - x^2))/((2*omega^2)/Cl^2 - omega^2/Cs^2 - ...
        2*x^2)))) + (1i*omega^2*(-2*lambda*x + 2*(lambda + 2*mu)*x + ...
     (omega^6*rhol)/(Cs^2*g*x*((2*omega^2)/Cl^2 - omega^2/Cs^2 - 2*x^2)*...
       sqrt(omega^2/Cl^2 - omega^2/Cp^2 - x^2)) - ...
     (4*omega^6*rhol*sqrt(omega^2/Cl^2 - omega^2/Cp^2 - x^2))/...
      (Cs^2*g*x*((2*omega^2)/Cl^2 - omega^2/Cs^2 - 2*x^2)^2) + ...
     (2*omega^6*rhol*sqrt(omega^2/Cl^2 - omega^2/Cp^2 - x^2))/...
      (Cs^2*g*x^3*((2*omega^2)/Cl^2 - omega^2/Cs^2 - 2*x^2)) - ...
     (4*mu*x*(omega^2/Cl^2 - x^2)*sqrt(omega^2/Cl^2 - omega^2/Cp^2 - x^2))/...
      (((2*omega^2)/Cl^2 - omega^2/Cs^2 - 2*x^2)*sqrt(omega^2/Cl^2 - omega^2/Cs^2 - ...
         x^2)) - (4*mu*x*(omega^2/Cl^2 - x^2)*sqrt(omega^2/Cl^2 - omega^2/Cs^2 - x^2))/...
      (((2*omega^2)/Cl^2 - omega^2/Cs^2 - 2*x^2)*sqrt(omega^2/Cl^2 - omega^2/Cp^2 - ...
         x^2)) - (8*mu*x*sqrt(omega^2/Cl^2 - omega^2/Cp^2 - x^2)*...
       sqrt(omega^2/Cl^2 - omega^2/Cs^2 - x^2))/((2*omega^2)/Cl^2 - omega^2/Cs^2 - ...
       2*x^2) + (16*mu*x*(omega^2/Cl^2 - x^2)*sqrt(omega^2/Cl^2 - omega^2/Cp^2 - x^2)*...
       sqrt(omega^2/Cl^2 - omega^2/Cs^2 - x^2))/((2*omega^2)/Cl^2 - omega^2/Cs^2 - ...
        2*x^2)^2)*((omega^2*rhol*sqrt(omega^2/Cl^2 - omega^2/Cp^2 - x^2))/...
      (Cs^2*((2*omega^2)/Cl^2 - omega^2/Cs^2 - 2*x^2)) + ...
     (1/g)*(lambda*(omega^2/Cl^2 - x^2) - (lambda + 2*mu)*(omega^2/Cl^2 - omega^2/Cp^2 - ...
         x^2) + (4*mu*(omega^2/Cl^2 - x^2)*sqrt(omega^2/Cl^2 - omega^2/Cp^2 - x^2)*...
         sqrt(omega^2/Cl^2 - omega^2/Cs^2 - x^2))/((2*omega^2)/Cl^2 - omega^2/Cs^2 - ...
         2*x^2))))/(x*(lambda*(omega^2/Cl^2 - x^2) - ...
      (omega^6*rhol*sqrt(omega^2/Cl^2 - omega^2/Cp^2 - x^2))/...
       (Cs^2*g*x^2*((2*omega^2)/Cl^2 - omega^2/Cs^2 - 2*x^2)) - ...
      (lambda + 2*mu)*(omega^2/Cl^2 - omega^2/Cp^2 - x^2) + ...
      (4*mu*(omega^2/Cl^2 - x^2)*sqrt(omega^2/Cl^2 - omega^2/Cp^2 - x^2)*...
        sqrt(omega^2/Cl^2 - omega^2/Cs^2 - x^2))/((2*omega^2)/Cl^2 - omega^2/Cs^2 - ...
        2*x^2))^2) + (1i*omega^2*((omega^2*rhol*sqrt(omega^2/Cl^2 - omega^2/Cp^2 - x^2))/...
      (Cs^2*((2*omega^2)/Cl^2 - omega^2/Cs^2 - 2*x^2)) + ...
     (1/g)*(lambda*(omega^2/Cl^2 - x^2) - (lambda + 2*mu)*(omega^2/Cl^2 - omega^2/Cp^2 - ...
         x^2) + (4*mu*(omega^2/Cl^2 - x^2)*sqrt(omega^2/Cl^2 - omega^2/Cp^2 - x^2)*...
         sqrt(omega^2/Cl^2 - omega^2/Cs^2 - x^2))/((2*omega^2)/Cl^2 - omega^2/Cs^2 - ...
         2*x^2))))/(x^2*(lambda*(omega^2/Cl^2 - x^2) - ...
     (omega^6*rhol*sqrt(omega^2/Cl^2 - omega^2/Cp^2 - x^2))/...
      (Cs^2*g*x^2*((2*omega^2)/Cl^2 - omega^2/Cs^2 - 2*x^2)) - ...
     (lambda + 2*mu)*(omega^2/Cl^2 - omega^2/Cp^2 - x^2) + ...
     (4*mu*(omega^2/Cl^2 - x^2)*sqrt(omega^2/Cl^2 - omega^2/Cp^2 - x^2)*...
       sqrt(omega^2/Cl^2 - omega^2/Cs^2 - x^2))/((2*omega^2)/Cl^2 - omega^2/Cs^2 - 2*x^2))) - 1i*h*sec(h*x)^2;
 
knn0=x-ff0/ffp0;
errr=abs(x-knn0)/knn0;
x=knn0;
end
beta_ii(i)=x;
beta_i(i)=real(beta_ii(i));
if imag(beta_ii(i))~=0
    beta_i(i)=nan;
end
val(i)=fff(beta_i(i));
left(i)=fffLeft(beta_i(i));
right(i)=fffRight(beta_i(i));
k_i(i)=sqrt(omega^2/Cl^2-beta_i(i)^2);
q_i(i)=sqrt(k_i(i)^2-omega^2/Cp^2);
s_i(i)=sqrt(k_i(i)^2-omega^2/Cs^2);
c_i(i)=omega/k_i(i);

 end




   if plott==1
 
  
      aa=0:pi/2/(h*20000):beta_ii(end);
for i=1:length(aa);
    avalLeft(i)=imag(fffLeft(aa(i)));
    avalRight(i)=imag(fffRight(aa(i)));
end
%here all the rh=(2n-1)pi/2 are excluded for the plot
RH=aa*h;
for i=1:floor(RH(end)/(pi))+1;
    [iii(i),jjj(i)]=find(abs(RH-((2*i-1)*pi/2))==min(abs(RH-((2*i-1)*pi/2))));
end
%here rr is excluded for the plot
avalLeft(jjj)=nan;

  
  %rigid bottom
[b0Rig,biRig]=acoustic_disp_rel(h,f,n_mod);
kRig=sqrt(omega^2/Cl^2-biRig.^2);
wRigreal=biRig(imag(kRig)==0);
wRigimag=biRig(real(kRig)==0);
valr=-omega^2./(g*aa);



width=900;  % Width of figure for movie [pixels]
height=600;  % Height of figure of movie [pixels]
left=200;     % Left margin between figure and screen edge [pixels]
bottom=200;  % Bottom margin between figure and screen edge [pixels]

figure
set(gcf,'Position', [left bottom width height])
plot([aa(1) aa(end)],[0 0],'--k');
hold on
plot(aa,avalLeft,'r');
hold on
plot(aa,avalRight,'b');
hold on
p22=plot(aa,valr,'--b');
hold on
%scatter(aa,aval,'x');
hold on
scatter(beta_ii,tan(beta_ii*h))
%ylim([-1.1*max(abs(tan(beta_i*h))) 1.1*max(abs(tan(beta_ii*h)))])

maxx=abs(tan(beta_i*h));
maxx(end+1:end+length(biRig))=abs(tan(biRig*h));
yup=max(maxx);
hold on  
%ylim([-1.1*max(abs(tan(RRROOOTT(1:n_mod,1)*h))) 1.1*max(abs(tan(RRROOOTT(1:n_mod,1)*h)))])
ylim([-1.1*yup 1.1*yup])



xlim([0 aa(end)])
set(gca,'XTick',[pi/2/h:pi/h:(2*n_mod-1)*pi/2/h],'XTickLabel',{'\pi/2','3\pi/2','5\pi/2','7\pi/2','9\pi/2','11\pi/2','13\pi/2','15\pi/2','17\pi/2'})
xlabel('r h')
ylabel('LHS/RHS')
 title(['h=',num2str(h),'m; \rho_l = ',num2str(rhol),'kgm-3; \rho_s = ',num2str(rhos),'kgm-3; C_l=',num2str(Cl),'ms-1; C_p=',num2str(Cp),'ms-1; C_s=',num2str(Cs),'ms-1; f=',num2str(f), 'Hz'],'fontsize',8)
 
%xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})
   end

