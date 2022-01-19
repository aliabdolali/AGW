function [beta0,beta_segn] = DM_acoustic_disp_rel(h,f,n_mod)
omega=2*pi*f;
gg=9.81;

stampo=0;
if stampo ==1
    figure    
end

% modo zero
L = gg/2/pi*(2*pi/omega)^2*tanh(sinh(omega*(h/gg)^0.5));
beta0 = 2*pi/(L);

% modi acustici
% n=1:2;
effe = @(x)abs(omega^2*h/9.81/x+tan(x));
for ir=1:n_mod
    effe_aprx = @(x)abs(omega^2*h/9.81/x+1/((n(ir)-1/2)*pi-x));
    dteta=pi/10000;
    kh=(ir-1)*pi+[pi/2:dteta:3/2*pi-dteta];
    left=omega^2.*h/9.81./kh;
    right=-tan(kh);
%     [x1,fval1]=fminbnd(effe_aprx,(ir-1)*pi+pi/2,(ir-1)*pi+pi/2+1e-2,optimset('TolX',1e-22,'Display','off'));
    [x,fval]=fminbnd(effe,(ir-1)*pi+pi/2,(ir-1)*pi+pi);
    x1 = -omega^2.*h.*pi*(ir-1/2)./(9.81-omega^2.*h); % (soluzione approssimata per x che tende a pi/2)
%     if f < 0.01317
        radice(ir)=x;
%     else
%         radice(ir)=x1;
%     end
    if omega^2.*h/9.81/((ir-1)*pi+pi/2) > 5
%          radice(ir)=(ir-1)*pi+pi/2;
    radice(ir) = x1;
    end
    if stampo ==1
        semilogy(kh,left)
        hold on
        semilogy(kh,right,'r')
        effep=abs(omega^2*h/9.81./kh+tan(kh));
        effepp=abs(omega^2*h/9.81./kh+(1/(n(ir)-1/2)*pi-x));
%         plot(kh,effep,'k')
        %             plot(kh,effepp,'g')
        grid on
        %             plot(kh,left-right,'r')
        %             plot(x1,omega^2*h/9.81/x1,'mo')
%         plot(x,omega^2*h/9.81/x,'ro')
        semilogy(radice(ir),omega^2*h/9.81/radice(ir),'k.','markersize',15)
    end
end
% ylim([0 omega^2.*h/9.81./(pi/2)*10])
beta_segn = radice./h;
if stampo == 1
    title(['h=' num2str(h) 'm, f=' num2str(f) 'Hz'])
end