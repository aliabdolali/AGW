function LLL=calct2(x0,y0,b,stepL,T,maxthreshL,minthreshL,maxerror,kaa1,kaa2,omegas1,omegas2,presionn1,presionn2)
    % this fuction helps to identify the right T
    const=presionn1/presionn2*(sqrt(kaa1)/sqrt(kaa2))*  abs((sin(kaa2*b)/sin(kaa1*b))) * abs((sin(omegas2*T)/sin(omegas1*T)));
    LLL=[];
    for L=[minthreshL:stepL:maxthreshL]

            a1=fcs(sqrt(kaa1./(pi.*x0)) .*  (L+y0));
            a2=fcs(sqrt(kaa1./(pi.*x0)) .*  (L-y0));

            a3=a1/2;
            a4=a2/2;

            a5=real(a3)+real(a4)+imag(a3)+imag(a4);
            a6=-real(a3)-real(a4)+imag(a3)+imag(a4);

            At=sqrt(a5.^2+a6.^2);

            a1=fcs(sqrt(kaa2./(pi.*x0)) .*  (L+y0));
            a2=fcs(sqrt(kaa2./(pi.*x0)) .*  (L-y0));

            a3=a1/2;
            a4=a2/2;

            a5=real(a3)+real(a4)+imag(a3)+imag(a4);
            a6=-real(a3)-real(a4)+imag(a3)+imag(a4);

            At2=sqrt(a5.^2+a6.^2);


            %M1=presion(2)/(const*abs(At));
            M1=const;
            M2=(At/At2);
            
            for i=1:length(M2)
                
                err=abs(abs(M1)-abs(M2(i)))/abs(M1)*100;
                
                if err>maxerror
                    LLL=[LLL, 0];
                else
                    LLL=[LLL, 1];
                end
                
            end
            
    end




end