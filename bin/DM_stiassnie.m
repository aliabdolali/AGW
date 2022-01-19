function[ETA,PT,ETA0,Pb0]=DM_stiassnie(h,b,tau,eta0,x,t,n)
%function[ETA,PT,ETA0,Pb0,PTSurface]=stiassnie10(h,b,tau,eta0,x,t,n)

g=10;
%h=4000;
%b=15000;
c=1500;
%tau=10;
%x=700000;
%t=100;
xb=x/h;
bb=b/h;
taub=tau*sqrt(g/h);
th=t*sqrt(g/h);
tb=th-0.5*taub;
cb=c/sqrt(g*h);

Pb0=1024*g*eta0*((2^(7/4))*sin((sqrt(2)*taub*((tb/xb)-1)^0.5)/2)*sin(sqrt(2)*bb*((tb/xb)-1)^0.5)*cos((sqrt(2)*xb*((tb/xb)-1)^(3/2))-(pi/4)))/(sqrt(pi)*taub*sqrt(xb)*((tb/xb)-1)^0.75*(2^1.5*sqrt((tb/xb)-1)+sinh((2^1.5)*sqrt((tb/xb)-1))));
ETA0=eta0*((2^(7/4))*sin((sqrt(2)*taub*((tb/xb)-1)^0.5)/2)*sin(sqrt(2)*bb*((tb/xb)-1)^0.5)*cosh((sqrt(2)*((tb/xb)-1)^(1/2)))*cos((sqrt(2)*xb*((tb/xb)-1)^(3/2))-(pi/4)))/(sqrt(pi)*taub*sqrt(xb)*((tb/xb)-1)^0.75*(2^1.5*sqrt((tb/xb)-1)+sinh((2^1.5)*sqrt((tb/xb)-1))));

for i=1:n
    al=(2*i-1)*pi/2;
  Pb(i)=1024*g*eta0*(2^2.5*sqrt(tb)*cb^1.5*(1-(xb/cb/tb)^2)^0.25*sin((al*cb*taub/2)/(sqrt(1-(xb/cb/tb)^2)))*sin((al*bb*xb/cb/tb)/(sqrt(1-(xb/cb/tb)^2)))*cos((al*sqrt(cb^2*tb^2-xb^2))+pi/4))/(sqrt(pi)*al^1.5*taub*xb);  
  
  eta(i)=-1*(2^2.5*sqrt(tb)*(1-(xb/cb/tb)^2)^1.25*sin((al*cb*taub/2)/(sqrt(1-(xb/cb/tb)^2)))*sin((al*bb*xb/cb/tb)/(sqrt(1-(xb/cb/tb)^2)))*cos((al*sqrt(cb^2*tb^2-xb^2))+pi/4))/(sqrt(pi)*al^2.5*taub*xb*cb^0.5);  
  
 
%    [beta0,beta_seg] = acoustic_disp_rel(h,(2*i-1)*c/4/h,i);
%    U(i)=beta_seg(i);
% PbSurface(i)=cos(U(i)*h)*Pb(i);


end
PT=sum(Pb);
%PTSurface=sum(PbSurface);
ETA=sum(eta);