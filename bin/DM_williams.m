function [eta0,eta,pressure0,pressure]= DM_williams(h,b,L,tau,zeta0,x,y,z,t,N, t_a)
% Evaluates bottom pressure WITH the effect of gravity as per BW /
% Kadri. Takes time (x,y,z,t) as argument. Outputs pressure.
    %global c rho W0 N U h
    c=1500;
    rho=1000;
    T = tau / 2;
   g=10;
   %acoustic modes
   %   for i=1:N
   %     U(i)=(2*i-1)*pi/2/h;
   % end
   
   


for i=1:N;
    [beta0,beta_seg] = DM_acoustic_disp_rel(h,(2*i-1)*c/4/h,i);
    U(i)=beta_seg(i);
end



    W0=zeta0/tau;
    pressure0 = zeros(1,length(t));
    p = zeros(N,length(t_a));
    eta0=zeros(1,length(t));
    eta=zeros(N,length(t_a));
    

    %zero mode
    Om_0 = sqrt(2*g/h)*sqrt((sqrt(g*h)*t./x)-1);
    K0=Om_0/sqrt(g*h);
    A00 = DM_A0(b,h,L,x,y,t);
    mu0=K0;
        auxg1 = (W0/pi/g) .* 8* A00 .* sin(K0*b) .* (sin(Om_0*T).*cosh(mu0.*h) ./ (K0 .* ((2.*mu0*h) + sinh(2.*mu0*h))));
        auxg2 = sqrt( (2*pi) ./ ((x/g)*sqrt(h/g)*Om_0 )  );
        auxg3 = cos(K0.*x - Om_0.*t + (pi/4));
    eta0 = auxg1 .* auxg2 .* auxg3;

    
    
     
   check = (1-(x./(c.*t_a)).^2);
    if check < 0
       eta = 0;
    else
         for n = 1:N
            Kn = DM_BigK(h, x, t_a, n);
          Om_n = DM_BigOmega(h, x, t_a,  n);
             A = DM_An(b,h,L,x, y, t_a,n);
             G = DM_Gn(b,h,tau,x, t_a, n);
             u = U(n);
             f = DM_omega(h,n);
          auxs1 = (W0/pi/g) .* A .* G .* ((2.*u.*Om_n).*cos(u.*h) ./ (Kn .* ((2.*u*h) + sin(2.*u*h))));
          auxs2 = sqrt( (2*pi) ./ (( x * f.^2) ./ (c * ( (Om_n.^2) - (f.^2) ).^(3/2) ) )  );
          auxs3 = cos(Kn.*x - Om_n.*t_a - (pi/4));
        et(n,:) = auxs1 .* auxs2 .* auxs3;
         end
         if N==1
             eta=et;
         else
      eta = sum(et);
         end
    end


  

    %%
    %Pressure
     %zero mode
    Om_0 = sqrt(2*g/h)*sqrt((sqrt(g*h)*t./x)-1);
    K0=Om_0/sqrt(g*h);
    A00 = DM_A0(b,h,L,x,y,t);
    mu0=K0;
        auxg1 = (rho*W0/pi) .* 8* A00 .* sin(K0*b) .* (sin(Om_0*T).*cosh(mu0.*z) ./ (K0 .* ((2.*mu0*h) + sinh(2.*mu0*h))));
        auxg2 = sqrt( (2*pi) ./ ((x/g)*sqrt(h/g)*Om_0 )  );
        auxg3 = cos(K0.*x - Om_0.*t + (pi/4));
    pressure0 = auxg1 .* auxg2 .* auxg3;


    %acoustic modes

   check = (1-(x./(c.*t_a)).^2);
    if check < 0
       pressure = 0;
    else
         for n = 1:N
            Kn = DM_BigK(h, x, t_a, n);
          Om_n = DM_BigOmega(h, x, t_a,  n);
             A = DM_An(b,h,L,x, y, t_a,n);
             G = DM_Gn(b,h,tau,x, t_a, n);
             u = U(n);
             f = DM_omega(h,n);
          aux1 = ((rho*W0)/pi) .* A .* G .* ((2.*u.*Om_n).*cos(u.*z) ./ (Kn .* ((2.*u*h) + sin(2.*u*h))));
          aux2 = sqrt( (2*pi) ./ (( x * f.^2) ./ (c * ( (Om_n.^2) - (f.^2) ).^(3/2) ) )  );
          aux3 = cos(Kn.*x - Om_n.*t_a - (pi/4));
          
          
        p(n,:) = aux1 .* aux2 .* aux3;
         end
         if N==1
             pressure=p;
         else
      pressure = sum(p);
         end
    end
  end

    
  
    
