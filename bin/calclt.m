function TTT=calclt(x0,y0,b,stepL,stepT,maxthreshT,minthreshT,maxthreshL,minthreshL,maxerror,kaa1,kaa2,omegas1,omegas2,presionn1,presionn2)

    const=presionn1/presionn2*(sqrt(kaa1)/sqrt(kaa2));%*  abs((sin(kaa2*b)/sin(kaa1*b)));
    TT=[];
    cn=0;
    TTT=[];
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
            M1=const*(At2/At);

            T1=minthreshT; % Period where I start iteration, minimun threshold
            error=5;
            while error>maxerror
                % La elevacion es 1 metros
                %T2=(1/M1)*abs(sin(omegas(2)*T1));
                cn=cn+1;
                M2=sin(omegas1*T1)/sin(omegas2*T1);

                error=abs((M1-M2))/M1*100;

                if error>maxerror
                    T1=T1+stepT; % Step for T, AT
                    TT(cn)=0;
                else
                    TT(cn)=T1;
                     T1=T1+stepT; % Step for T, AT
                     error=5;
                     
                end
                if T1>maxthreshT % Period where i stop iterating, maximun threshold
                    break
                end
            end


            TTT=[TTT;TT];
            TT=[];
            cn=0;
    end




end