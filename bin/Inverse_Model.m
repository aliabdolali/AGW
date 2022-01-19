
function [b_m, L_m, T_m, W0_m, theta, X0, Y0] = Inverse_Model_2022_01_10(min_f, max_f, npts, nnn, Fs, h, r_dist, DT, Pmax, clas, magnitude, TT, YY, ii, alpha)

b_final = nan; L_final=nan; T_final=nan; W0_final=nan; theta=nan; X0=nan; Y0=nan;

%IM_input_signal_segment

% for ii=2:floor(TT(end)/DT)

figure(1);
%clf      
clear t
clear y
t(:,1)=TT(TT<=ii*DT);
y(:,1)=YY(TT<=ii*DT);
torig=t;

m=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:floor(t(end)/DT)
    clear ta1
    clear ya1
   ta1=t(t<=i*DT & t>=(i-1)*DT);
   ya1=y(t<=i*DT & t>=(i-1)*DT);
   if max(abs(ya1))>=Pmax
       tm(m+1,1)=ta1(1);
       tm(m+2,1)=ta1(end);
       m=m+2;
   end
end

if m>=1
tm=unique(tm); %to make them continues
clear tmm
tmm(1)=tm(1);
if length(tm)>=3
for iii=2:length(tm)-1
    if tm(iii)+DT~=tm(iii+1) & tm(iii)-DT~=tm(iii-1)
   tmm(end+1,1)=tm(iii);
    end
    tmm(end+1)=tm(end);       
end
else
tmm=tm; 
end
clear yy
clear tta
for j=1:length(tm)-1
    clear ytmp
    clear ttmp
    ytmp(:,1)=y(t<=tm(j+1) & t>=tm(j));
    ttmp(:,1)=t(t<=tm(j+1) & t>=tm(j));
    if j==1
        yy(:,1)=ytmp;
        tta(:,1)=ttmp;
    else
        yy(end+1:end+length(ytmp),1)=ytmp;
        tta(end+1:end+length(ytmp),1)=ttmp; 
    end
end
end


%% Spectral content
len=25000;
s = spectrogram(y,len,[],len,Fs);

%% Here I chose the maximum possible 1st mode frequency
%%disp('Select maximum observed potential 1st mode frequency in the spectrogram and press ENTER');
%wg = warndlg('Zoom in the spectrogram and click OK when ready, then select the maximum observed potential 1st mode frequency in the spectrogram and press ENTER', 'Zoom');
%wg.Position(1) = wg.Position(1)-00;
%wg.Position(2) = wg.Position(2)-250;
%waitfor(wg);
%[xf,yf] = ginput;
%max_f=yf;
%% Here I chose the minimum possible 1st mode frequency
%%disp('Select minimum observed potential 1st mode frequency in the spectrogram and press ENTER');
%wg = warndlg('Zoom in the spectrogram and click OK when ready, then select the minimum observed potential 1st mode frequency in the spectrogram and press ENTER', 'Zoom');
%wg.Position(1) = wg.Position(1)-00;
%wg.Position(2) = wg.Position(2)-250;
%waitfor(wg);
%[xf,yf] = ginput;
%min_f=yf;
%% RANGES FOR PROPERTIES BASED ON WELLS AND COPPERSMITH 1994 RELATED TO MW
% USING WELLS 1994 
T_range=(5.25: 0.5: 50.25); %[s]  THIS VALUE REMAINS FOR BOTH APPROACHES

MW=magnitude;
a=-3.22 ;b=0.69 ; %a=-3.55 ;b=0.74 ;
SRL=(10^(a+b*MW));
L_range=(round(SRL*0.5)*10^3: 1*10^3: round(SRL*1.75)*10^3); %[m]
% WITH THE RANGES WE REPRESENT THE SCATTER FOUND ON THE PARAMETERIZATIONS

a=-3.49 ;b=0.91 ; 
RA=10^(a+b*MW);

BRL=RA/SRL;
b_range=(round(BRL*0.5)*10^3: 1*10^3: round(BRL*1.75)*10^3); %[m]

% a=-5.46 ;b=0.82 ;
a=-4.80; b=0.69 ;
DISP=10^(a+b*MW);

L_range=L_range./2; b_range=b_range./2; % We work with half length and width
stepL=1000; stepT=0.2; stepb=1*10^3; % Setup parameter! 

%% Set up parameters for the model
% I assume I know the frequency range after looking at the spectrogram
d_step=100*10^3; % Grid step for considered X0 solutions
frec=Fs; c=1500; ro=1000;
window_ener=4000; % Size of the short time energy analysis window in samples
n_iterations_pks=3; % Number of iterations for envelope tracking
b_sol_point=20; % Number of solutions asked for sin(kb)
maxerror=0.3; % maximun error that can have the solution, this is for T and L
percent_diff_T_sol=20; % Percentage to consider solutions from histogram for T
%% HERE I CHOSE THE START AND END OF MY SIGNAL
% Short time energy is computed
ww=y; power=[];
for i =1:length(t)/window_ener
    aux2= y(1 +window_ener*(i-1) : window_ener+window_ener*(i-1)).* y(1 +window_ener*(i-1) : window_ener+window_ener*(i-1));
    power=[power, sum(aux2)];
   tt(i)=t((i-1)*window_ener+1);
end

    figure(1);
    hold off
    subplot(4,2,1);
    plot(t,ww,'b')
    %title('IMPORTANT TO CLICK ON THE SIGNAL CHART AND NOT ON THE ENERGY DISTRIBUTION CHART!')
    xlabel('Time (s)')
    ylabel('Pressure (Pa)')
    axis([min(t)+30 max(t)-30 min(ww) max(ww)])
    ylim([-60 60])
    hold on
    grid on
    subplot(4,2,2);

  
[~,F,T,P] =spectrogram(y,len,[],len,Fs,'yaxis');
imagesc(T, F, 10*log10(P+eps)); % add eps like pspectrogram does
axis xy
ylabel('Frequency (Hz)')
xlabel('Time (s)')
colorbar;
h1.Label.String = 'Power/frequency (dB/Hz)';
h1.Location='northoutside';

if m>=1
hold on
for jj=1:length(tmm)
    plot([tmm(jj) tmm(jj)],[min(F(:)) max(F(:))],'--k')
    hold on
end
end
xlim([t(1)+30,t(end)-30])
    subplot(4,2,3);
    if m>=1
    plot(tt,power,'b')
    xlabel('Time (s)')
    ylabel('Power [RMS]')
    axis([torig(1)+30 torig(end)-30 min(power) max(power)])
 
    else
      plot(torig,nan*ones(size(torig)),'b')
      xlim([torig(1) torig(end)])
          xlabel('Time (s)')
    ylabel('Power [RMS]')
    end
            box on
         axis on
         
         
        subplot(4,2,4);
         plot(torig,nan*ones(size(torig)),'b')
         xlim([torig(1)+30 torig(end)-30])
        box on
         axis on
         xlabel('Time (s)')
ylabel('Frequency [Hz]')
title('Potential frequency distributions for the first acoustic mode')

        
%     x0p=300;  y0p=200;
%     widthp=1000;  heightp=1000;
%     set(gcf,'position',[x0p,y0p,widthp,heightp])

if m>=1
        
        %close(figure(2))%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    buttonSelections = 1; count2=0;
    
    while buttonSelections == 1
        
        
        
        
        %% Here a FIRST point in the plot is defining the BEGGINING of the region-----
        %%disp('Click OK. Then, select the beggining of the disturbance in the signal and press ENTER');
        %wg = warndlg('Zoom in the signal and click OK when ready. Then, select the beggining of the disturbance in the signal and press ENTER', 'Zoom');
        %wg.Position(1) = wg.Position(1)-00;
        %wg.Position(2) = wg.Position(2)-250;
        %waitfor(wg);
        
        %[x,y] = ginput;
        
        x1=tta(1);
        x=x1;
        aa=[];
        for i=1:length(x)    
            aa(i,:)=abs(x(i)-t); 
        end
        locc=[];
        for i=1:length(x)     
            locc(i)=find(aa(i,:)==min(aa(i,:)));
        end
        
        %% Here a second point in the plot is defining the end of the region-----
        %%disp('Click OK. Then, select the endof the disturbance in the signal and press ENTER');
        %%wg = warndlg('Zoom in the signal and click OK when ready. Then, select the end of the disturbance in the signal and press ENTER', 'Zoom');
        %%wg.Position(1) = wg.Position(1)-300;
        %wg.Position(1) = wg.Position(1)-00;
        %wg.Position(2) = wg.Position(2)-250;
        %waitfor(wg);
     
        %[x,y] = ginput;
        x2=tta(end);
        x=x2;
        aa=[];
        for i=1:length(x)    
            aa(i,:)=abs(x(i)-t); 
        end
        locc2=[];
        for i=1:length(x)     
            locc2(i)=find(aa(i,:)==min(aa(i,:)));
        end
        
        %% THIS NEED A CHECK CAUSE IM NOT SURE IF IT WORKS FOR MORE THAN 1 AREA
       % promptMessage = sprintf('Do you want more areas');
       % titleBarCaption = 'Yes or No';
       % button = questdlg(promptMessage, titleBarCaption, 'Yes', 'No', 'Yes');
        
       %   if strcmpi(button, 'Yes')
       %     buttonSelections= 1;
       %   else
            buttonSelections = -1;
       %   end
    end
    
    
    
    t_length=t(locc2)-t(locc);% Approximated length in seconds of the audio content
    Audio_start=t(locc);% Beggining of audio content relative to the recorded signal
    
    
    
    %% NEW IMPLEMENTATION LOCATION
    
    signal_loc=ww(locc:locc2);
    signal_t=t(locc:locc2);
    
    pks_aux=[]; locs_aux=[];locs_aux2=[];
    pks=signal_loc;
    [pks,locs_aux]=findpeaks(pks);
    locs_aux2=locs_aux;
    for itera=1:n_iterations_pks % n iterations +1 which is out of the loop
    
        [pks,locs_aux]=findpeaks(pks);
        locs_aux2=locs_aux2(locs_aux);
    end
    tims_new=signal_t(locs_aux2); pks_new=pks;
    
    periods=tims_new(2:end)-tims_new(1:end-1);
    frequencies=(2*pi)./periods;
    
    % check that the frequency is decaying and apply formulation
    final_freq=[];
    period_aux=(tims_new(2)-tims_new(1));
    periods=[]; tims_new2=[]; pks_new2=[];
    periods=[periods, period_aux];
    tims_new2=[tims_new2, tims_new(1)];
    pks_new2=[pks_new2, pks_new(1)];
    counter=0;
    for i=2:(length(pks_new)-1)
        
        ppp=(tims_new(i+1)-tims_new(i-counter));
        
        if ppp>period_aux
            periods=[periods, ppp];
            tims_new2=[tims_new2, tims_new(i)];
            pks_new2=[pks_new2, pks_new(i)];
            
            period_aux=ppp;
            counter=0;
        else
            period_aux=period_aux;
            counter=counter+1;
        end
        
    end
    
    
    omegas_new=(2.*pi)./periods;
    X0=[]; Y0=[]; T0=[]; 
    w=pi*c/(2*h);
    % periods; tims_new2; pks_new2;
    for iii=1:length(tims_new2)
        nlocs3(iii,:)=sort(randperm(length(tims_new2),2));
        nn1=nlocs3(iii,1);
        nn2=nlocs3(iii,2);
        % LOCATION AND ERUPTION TIME
        % X0, Y0 and t0 are retrieved and saved in this piece of code
        aux1 = (1 - ( (pi*c)/(2*h*omegas_new(nn1)) )^2 )^(-1/2);
        aux2 = (1 - ( (pi*c)/(2*h*omegas_new(nn2)) )^2 )^(-1/2);
        x02= (tims_new2(nn2)-tims_new2(nn1))*c/ (aux2-aux1); % retrieved X coordinate of the hydrophone
    
        aux1=( 1 -  (pi*c/(2*h*omegas_new(nlocs3(1))))^2 )^(-1/2);
        t02= (tims_new2(nn1)) - (x02/c) * aux1; % retrieved eruption time
    
        yy=abs(sqrt(abs(t02*c)^2-x02^2));
        x0=x02;
        y0=yy;
        
        if w<omegas_new(nn1) & w<omegas_new(nn2)
            X0=[X0, x0];
            Y0=[Y0, y0];
            T0=[T0, t02];
        end
        
    end
    X0_new=mean(X0);
    Y0_new=mean(Y0);
    tot_dist=sqrt(X0_new^2+Y0_new^2);
    T0_new=mean(T0);
    
    
    
    
    %%
    
    
    
    %% THE POTENTIAL FREQUENCY DISTRIBUTIONS RELATED TO THE POTENTIAL LOCATIONS ARE COMPUTED
    t_r=t; w=pi*c/(2*h); % Frequency for the FIRST MODE
    x00=0:d_step:r_dist; y00=sqrt((r_dist).^2-(x00).^2);
    t0=0-(sqrt(y00(1).^2+x00(1).^2) ./ c);
    OM=[];
    for i=1:length(x00)
    
        t=[-t0  -t0+t_length]; %%%% THIS IS THE LENGTH OF SIGNAL WE WANT TO CONSIDER
        omega=(w ./ sqrt(1- (x00(i)./(c.*t)).^2 ) ); % Frequency
        OM=[OM; omega];
    end
    
    pos=[];
    % Verify that the frequencies are in the required range
    for i=1:length(x00)
        if OM(i, 1)<max_f && OM(i, 2)>min_f
            pos=[pos, i];
        end
    end
    
    
    X0_p=x00(pos); Y0_p=y00(pos);
    omega=[];
    for il=1:length(pos)  
        x0=x00(pos(il));
        y0=y00(pos(il));
        t=[-t0 : 1/frec : -t0+t_length];
        omega(il,:)=(w ./ sqrt(1- (x0./(c.*t)).^2 ) );
        hold on
        %plot(t-(t(1))+Audio_start,omega(il,:),'k','Linewidth',1)
    end
    
    
    t1=t_r(locc:locc2); y1=ww(locc:locc2);
    %% Tracking the envelope
    
    pks_aux=[]; locs_aux=[];locs_aux2=[];
    pks=y1;
    [pks,locs_aux]=findpeaks(pks);
    locs_aux2=locs_aux;
    for itera=1:n_iterations_pks % n iterations +1 which is out of the loop
    
        [pks,locs_aux]=findpeaks(pks);
        locs_aux2=locs_aux2(locs_aux);
    end
    tims2=t1(locs_aux2); pks2=pks;
    
    periods=tims2(2:end)-tims2(1:end-1);
    frequencies=(2*pi)./periods;
    
    
    
    time1=(1:length(omega))/ frec; % times associated to omegas
    time2=tims2-t_r(locc);% times associated to pks
    finaltims=time2; finalpks=pks2;
    time22=tims2;
    % find the location of the related omegas to the tracked envelop
    aa=[];
    for i=1:length(time2)    
        aa(i,:)=abs(time2(i)-time1); 
    end
    locc=[];
    for i=1:length(time2)     
        locc=[locc, find(aa(i,:)==min(aa(i,:)))];
    end
    omeg=[];
    for ili=1:size(omega,1)
        omeg(ili,:)=omega(ili,locc);
    end
    
    
    % These are the points that the model will use
    figure(1);
    subplot(4,2,1);
    hold on
    plot(t1,y1,'r')
    hold on 
    plot(time22, pks2, '.m')
    xlabel('Time (s)')
    ylabel('Pressure [Pa]')
    title('Points used by the model')
    %xlim([min(time2), max(time2)])
    
    subplot(4,2,4);
    hold on
    for ili=1:size(omega,1)
        plot(time1+x1, omega(ili,:),'k')
        hold on
        plot(time22, omeg(ili,:),'ro')
    end
    xlabel('Time (s)')
    ylabel('Frequency [Hz]')
    title('Potential frequency distributions for the first acoustic mode')
    %xlim([min(time2), max(time2)])
    xlim([torig(1)+30, torig(end)-30])
    box on
    axis on
    
    
    
    % x0p=100;  y0p=220;
    % widthp=1000;  heightp=1000;
    % set(gcf,'position',[x0p,y0p,widthp,heightp])
    
      
    %% NOW ORIENTATION, LOCATION AND FREQUENCY DISTRIBUTION ARE KNOWN
    %prompt = 'How many solutions [1,Infty) per potential first mode frequency distribution do you want? (optimal 5-20) ';
    %nnn = input(prompt) % NUMBER OF SETS OF SOLUTIONS THAT I WANT
    
    %nnn=10;% Here nnn number of non repeated combinations from the available points
    %% are made they are also sorted getting the latest points at the beggining
    %% since those are the selected to compute the location and eruption time
    
    %prompt = 'How many points [4,max_number_of_points] per solution should be used by the model? (optimal 4-6)';
    %npts = input(prompt) % NUMBER OF POINTS FROM THE SIGNAL THAT I WANT FOR COMPUTING EACH SET OF SOLUTIONS
    %nlocs=[];
    %npts=4;
    %% I obtain solutions for all possible omega distributions
    ome_ga=omega;
    W0S=[]; TS=[]; LS=[]; BS=[]; W0S1=[];
    tic
    for jj=1:size(ome_ga,1)
        
        x0=X0_p(jj); y0=Y0_p(jj);
        omegs=omeg(jj,:);
        for i=1:nnn
           %This part ensures that the points are not too close
           nlocs(i,:)=sort(randperm(length(finaltims),npts)) ;
    
           tiempo1=finaltims(nlocs(i,1));
           tiempo2=finaltims(nlocs(i,2));
           while abs(abs(tiempo1)-abs(tiempo2))<10 % TUNING PARAMETER IMPORTANT DISTANCE IN SECONDS BUT OPNLY FROM THE 2 FIRST POINTS, THEY COMPUTE LOCATION
                nlocs(i,:)=sort(randperm(length(finaltims),npts)) ;
                tiempo1=finaltims(nlocs(i,1));
                tiempo2=finaltims(nlocs(i,2));
           end
        end
        %% Properties of the slender fault
        X00=[]; Y00=[]; T00=[];q=0; bsol=[]; Tsol=[]; W00sol=[]; Lsol=[]; 
        pruebaT=[]; iteration=0; W0prue=[]; W0prue2=[]; 
        t02=sqrt(x0^2+y0^2)/c;% t0 is defined because the location is known
        for j=1:length(nlocs)
            % Points from and related to the signal are defined
            iteration=iteration+1;
            tiempos=finaltims(nlocs(j,:))+t02;
            presiones=finalpks(nlocs(j,:));  
            omegas=omegs(nlocs(j,:));
            kaa=omegas.*x0./(c^2.*(tiempos+t02)); % Wavenumber
    
            %% b (WIDTH)
    
            n=[1:b_sol_point]; % I want the first X solutions for each point in the analytical solution of 1=sin(kb)
            b2=[];
            for i = 1:length(tiempos)
                b2(i,:)= ((n-1/2)*pi)/kaa(i);
            end
            bb=b_range; histo=[];
    
            % We track around what points is the convergence in the solutions
            for j= 1:length(tiempos)
                for i=1:length(bb)-1
                    range=bb(i)+stepb;
                    auxil=find(b2(j,:)<range & b2(j,:)>bb(i));
                    rr=isempty(auxil);    
                    if j==1
                        if rr==1
                            histo=[histo,0];
                        else
                            histo=[histo,1];
                        end
                    else
                        if rr==1
                            histo(i)=histo(i);
                        else
                            histo(i)=histo(i)+1;
                        end
                    end
                end
            end
    
            % The convergence of points is found
            B=find(histo==max(histo));
            b=bb(B)+0.5*10^3;
            % If more than a solution is found they are averaged
            if length(b)>1
                baux=sum(b)/length(b);
                b=baux;
            end
    
            %% HERE T IS CALCULATED EVALUATION ALL POSSIBLE SOLUTIONS FOR EACH L IN A RANGE THAT WE ESTABLISH
            maxthreshT=max(T_range); %Range of possible T
            minthreshT=min(T_range);
            maxthreshL=max(L_range); %Range of possible L
            minthreshL=min(L_range); 
    
            % Combinations between the chosen points in the chart because we need 2 points to compute the formulation
            combos = nchoosek(1:length(tiempos),2);
            L=[minthreshL:stepL:maxthreshL];
            TTT=[];
            % All possible T are evaluated
            for i=1:size(combos,1)
                n1=combos(i,1); n2=combos(i,2);
                kaa1=kaa(n1); kaa2=kaa(n2);
                omegas1=omegas(n1); omegas2=omegas(n2);
                presionn1=presiones(n1); presionn2=presiones(n2);
                TTT(:,:,i)=calclt(x0,y0,b,stepL,stepT,maxthreshT,minthreshT,maxthreshL,minthreshL,maxerror,kaa1,kaa2,omegas1,omegas2,presionn1,presionn2);
            end
    
            Ts=(minthreshT:stepT:maxthreshT);
            histo=[]; xx=[]; rr=[];
            % The convergence of solutions is evaluated
            for i=1:size(combos,1)
                for j=1:length(Ts)
                    xx=find(TTT(:,:,i)>Ts(j) & TTT(:,:,i)<Ts(j)+stepT);
                    rr=isempty(xx);
                    histo(i,j)=length(xx);   
                end
            end
            fff=sum(histo);
    
            % THE UNCERTAINTY OF T IS EVALUATED AND ASSOCIATED TO THIS SOLUTION
            sum2=fff;
            val=[]; ind=[];
            [val ind] = sort(sum2,'descend');
    
            % I ONLY CONSIDER THE SOLUTIONS WITH LESS THAN percent_diff_T_sol DIFFERENCE TO
            % THE HIGHEST CONVERGENCE OF SOLUTIONS
            top2=[val(1)];
            TTS=[Ts(ind(1))]+stepT/2;
            if length(val)>1
                if (val(1)-val(2))/val(1)*100<percent_diff_T_sol 
                    top2=[val(1),val(2)];
                    TTS=[Ts(ind(1)), Ts(ind(2))]+stepT/2;
                end
                if (val(1)-val(3))/val(1)*100<percent_diff_T_sol
                    top2=[val(1),val(2),val(3)];
                    TTS=[Ts(ind(1)), Ts(ind(2)), Ts(ind(3))]+stepT/2;
                end
                if (val(1)-val(4))/val(1)*100<percent_diff_T_sol
                    top2=[val(1),val(2),val(3),val(4)];
                    TTS=[Ts(ind(1)), Ts(ind(2)), Ts(ind(3)), Ts(ind(4))]+stepT/2;
                end
            end
    
            pruebaT=[pruebaT, mean(TTS)];
            % Here the chosen possible solutions for T are evaluated again by
            % looking at the error that the generate using the function calct2
            cnr=[]; T=[];
            for i=1:size(combos,1)
                n1=combos(i,1); n2=combos(i,2);
                kaa1=kaa(n1); kaa2=kaa(n2);
                omegas1=omegas(n1); omegas2=omegas(n2);
                presionn1=presiones(n1); presionn2=presiones(n2);
    
                for j=1:length(TTS)
                   LLL=calct2(x0,y0,b,stepL,TTS(j),maxthreshL,minthreshL,maxerror,kaa1,kaa2,omegas1,omegas2,presionn1,presionn2); 
                   auxill=find(LLL==1);
                   cnr(j)=length(auxill);
                end
                indi=find(cnr==max(cnr));
                T=[T, TTS(indi)];
            end
    
            cnr=[];
            for i=1:length(TTS)
                tb=TTS(i);
                tb1=find(T==tb);
                cnr=[cnr, length(tb1)];
            end
            T=TTS(find(cnr==max(cnr)));
    
            if length(T)>1
                Taux=sum(T)/length(T);
                T=Taux;
            end
            Tr=T; % The convergence is found around certain period
    
            %% W0
            dmax=DISP*2;
            if dmax>9
                dmax=9; % Observes in Wells (1994) dataset
            end
            W_range=(round((DISP*0.5)/(Tr*2),3): 0.01: round((dmax)/(Tr*2),3)); %[m/s]
    
            LL=L_range; W00=W_range; T=Tr;
            % redefine k1, omega and t, we just need 2 points n1 and n2
            n1=2; n2=3;
    
            k1=[kaa(n1), kaa(n2)];
            omega=[omegas(n1), omegas(n2)];
            t=[tiempos(n1), tiempos(n2)];
            pressure=[presiones(n1), presiones(n2)];
            solucion1=[]; solucion2=[]; sol1=[];  sol2=[];
    
            % Two surfaces are generated which display the error for each combination of L and W0 when 
            % they generate pressure compared to the measured pressure, since the
            % solution of W0 is unique there is a line for similar values of W0
            % that is averaged in order to give a final solution for W0
    
            for i=1:length(LL)
                for j=1:length(W00)
    
                    E=b/LL(i);
                    l=E*LL(i); % b
                    X=x0*(E^2); Y=y0*E;
                    Y1pos= (l+Y)/2 ; Y1neg= (l-Y)/2 ;
                    v=X./k1;
                    X1=v./2;    
                    a1=fcs(sqrt(2./(pi.*X1)) .*  Y1pos);
                    a2=fcs(sqrt(2./(pi.*X1)) .*  Y1neg);
                    a3=a1/2; a4=a2/2;
                    a5=real(a3)+real(a4)+imag(a3)+imag(a4);
                    a6=-real(a3)-real(a4)+imag(a3)+imag(a4);
                    At=sqrt(a5.^2+a6.^2); % Envelope
    
                    p1= ro.*W00(j).*At  .*   ((2.^(5./2) .* (c.^3) .*  ((t).^(1/2)))  ./  (h.* pi.^(1./2) .* w.^(3./2) .*x0 ) );
                    p2= (1 - (x0./ (c.*(t))).^2).^(1./4);
                    p3= sin( k1.* b  );
                    p4= sin ( omega .* T   );
                    aux= p1.*p2.*p3.*p4;
                    press=abs(aux); % Pressure signal       
    
                    val1=pressure(1)-press(1);
                    val2=pressure(2)-press(2);
                    solucion1=[solucion1, abs(val1)];
                    solucion2=[solucion2, abs(val2)];
                end
                sol1=[sol1; solucion1];
                sol2=[sol2; solucion2];
                solucion1=[]; solucion2=[];
            end
    
            aa1=sum(sol1); aa2=sum(sol2);
            ss1=find(aa1==min(aa1));
            ss2=find(aa2==min(aa2));
            W03=(W00(ss1)+W00(ss2))/2;
            
            %%%% FIND THE MINIMUM OF EACH SURFACE
            [row,col1] = find(sol1==min(sol1(:)));
            [row,col2] = find(sol2==min(sol2(:)));
            W04=(W00(col1)+W00(col2))/2;
            %% L
    
            LLL=[];
            for i=1:length(kaa)
                LL=L_range;
    
                E=b./LL;
                l=E.*LL; % b
                X=x0.*(E.^2); Y=y0.*E;
                Y1pos= (l+Y)./2 ; Y1neg= (l-Y)./2 ;
                v=X./kaa(i);
                X1=v./2;    
                a1=fcs(sqrt(2./(pi.*X1)) .*  Y1pos);
                a2=fcs(sqrt(2./(pi.*X1)) .*  Y1neg);
                a3=a1/2; a4=a2/2;
                a5=real(a3)+real(a4)+imag(a3)+imag(a4);
                a6=-real(a3)-real(a4)+imag(a3)+imag(a4);
                At1=sqrt(a5.^2+a6.^2);
    
                p11= ro.*W04.*At1  .*   ((2.^(5./2) .* (c.^3) .*  ((tiempos(i)).^(1/2)))  ./  (h.* pi.^(1./2) .* w.^(3./2) .*x0 ) );
                p22= (1 - (x0./ (c.*(tiempos(i)))).^2).^(1./4);
                p33= sin( kaa(i).* b  );
                p44= sin ( omegas(i) .* T   );
    
                Ataux=p11.*p22.*p33.*p44;
    
                min1=(abs(abs(presiones(i))-abs(Ataux(:))));
                loc1=find(min1==min(min1));
                LLL=[LLL, LL(loc1)]; %FIND THE L WITH THE SMALLEST ERROR FOR EVERY K
            end
    
            L=sum(LLL)/length(LLL);
    
            %% W0
            frec=1; % Frequency for the signals that i generate to calculate W0
            W00=W_range;
    
            t=[r_dist/c+1 : 1/frec : t_length+r_dist/c+1];
            maxpres=[];
            for jj =1:length(W00)
                k1=(w.*x0)  ./  ((c.^2) .* (t) .* sqrt(1-   ((x0./  (c.*(t))).^2))); % Wave number
                omega=(w ./ sqrt(1- (x0./(c.*t)).^2 ) ); % Frequency
    
                E=b/L;
                l=E*L; % b
                X=x0*(E^2); Y=y0*E;
                Y1pos= (l+Y)/2 ; Y1neg= (l-Y)/2 ;
                v=X./k1;
                X1=v./2;    
                a1=fcs(sqrt(2./(pi.*X1)) .*  Y1pos); 
                a2=fcs(sqrt(2./(pi.*X1)) .*  Y1neg);
                a3=a1/2; a4=a2/2;
                a5=real(a3)+real(a4)+imag(a3)+imag(a4);
                a6=-real(a3)-real(a4)+imag(a3)+imag(a4);
                At=sqrt(a5.^2+a6.^2); % Envelope
    
                p1= ro.*W00(jj).*At  .*   ((2.^(5./2) .* (c.^3) .*  ((t).^(1/2)))  ./  (h.* pi.^(1./2) .* w.^(3./2) .*x0 ) );
                p2= (1 - (x0./ (c.*(t))).^2).^(1./4);
                p3= sin( k1.* b  );
                p4= sin ( omega .* T   );
                aux= p1.*p2.*p3.*p4;
    
                maxpres=[maxpres, max(abs(aux))];
            end
    
            difer=abs(maxpres-max(pks2));
            difer2=W00(find(difer==min(abs(maxpres-max(pks2)))));
            if length(difer2)==1
                W02=[ W00(find(difer==min(abs(maxpres-max(pks2)))))];
            else
                W02=[ W04];
            end
    
            bsol=[bsol, b]; Tsol=[Tsol, T]; Lsol=[Lsol, L];
            W00sol=[W00sol, W02];
            %W0prue=[W0prue, W03];
            W0prue2=[W0prue2, W04];
    
        end
        
        W0S1=[W0S1; W00sol];
        W0S=[W0S; W0prue2];
        TS=[TS; Tsol]; LS=[LS; Lsol]; BS=[BS; bsol];
        
    end
    toc
    % sum(sum(W0S1))/(size(W0S1,1)*size(W0S1,2))
    % sum(sum(W0S))/(size(W0S,1)*size(W0S,2))
    % W0S=W0S1;
    
    %%  CHARTS FOR THE SOLUTIONS
    
    LH=L_range; TH=T_range; bH=b_range;
    
    W0_final=[];
    for i=1:size(ome_ga,1)    
        W0_final=[W0_final, W0S(i,:)];
    end
        
    T_final=[];
    for i=1:size(ome_ga,1)    
        T_final=[T_final, TS(i,:)];
    end
    
    L_final=[];
    for i=1:size(ome_ga,1)    
        L_final=[L_final, LS(i,:)];
    end
    
    B_final=[];
    for i=1:size(ome_ga,1)    
        B_final=[B_final, BS(i,:)];
    end
    
    b_final=B_final; 
    LH=L_range; TH=T_range; bH=b_range;
    xbins=(min(TH)-stepT:stepT:max(TH)+stepT);
    
    histo=[];
    for i=1:length(xbins)-1
        
        auxil=find(T_final(:)<=xbins(i+1) & T_final(:)>=xbins(i));
        rr=isempty(auxil) ;
        histo(i)=length(auxil);
    end
    
    bins=xbins(1:(end-1))+(stepT/2);
    xi = linspace(min(bins), max(bins), 350);         
    yi = interp1(bins, histo, xi, 'v5cubic');
    hhh1=figure(1);
    hold off
    hb8=subplot(4,2,8);
    hold off
    a=area(xi, yi);
    a.FaceAlpha = 0.25;
    yl1=[0 30];
    vl=[mean(T_final) mean(T_final)];
    hold on
    plot(vl,yl1,'r--')
    xlim([min(bins)-stepT, max(bins)+stepT])
    ylim([0,max(histo)+2])
    xlabel('T [s]')
    ylabel('# solutions')
    title(['Average = ',num2str(vl(1)), ' s'])
    
    
    %%%%%%%%%%%%%%%%%%%
    
    xbins=(min(LH)-stepL:stepL:max(LH)+stepL);
    histo=[];
    for i=1:length(xbins)-1
        
        auxil=find(L_final(:)<=xbins(i+1) & L_final(:)>=xbins(i));
        rr=isempty(auxil) ;
        histo(i)=length(auxil);
    end
    
    bins=xbins(1:(end-1))+(stepL/2);
    xi = linspace(min(bins), max(bins), 350);         
    yi = interp1(bins, histo, xi, 'v5cubic');
    subplot(4,2,7);
    hold off
    a=area(xi, yi);
    a.FaceAlpha = 0.25;
    
    yl1=[0 30];
    vl=[mean(L_final) mean(L_final)];
    hold on
    plot(vl,yl1,'r--')
    
    xlim([min(bins)-stepL, max(bins)+stepL])
    ylim([0,max(histo)+2])
    xlabel('L [m]')
    ylabel('# solutions')
    title(['Average = ',num2str(vl(1)), ' m'])
    
    %%%%%%%%%%%%%%%%%%%
    
    xbins=(min(bH)-stepb:stepb:max(bH)+stepb);
    histo=[];
    for i=1:length(xbins)-1
        auxil=find(b_final(:)<=xbins(i+1) & b_final(:)>=xbins(i));
        rr=isempty(auxil) ;
        histo(i)=length(auxil);
    end
    
    bins=xbins(1:(end-1))+(stepb/2);
    xi = linspace(min(bins), max(bins), 350);         
    yi = interp1(bins, histo, xi, 'v5cubic');
    subplot(4,2,5);
    hold off
    a=area(xi, yi);
    a.FaceAlpha = 0.25;
    
    yl1=[0 30];
    vl=[mean(b_final) mean(b_final)];
    hold on
    plot(vl,yl1,'r--')
    
    xlim([min(bins)-stepb, max(bins)+stepb])
    ylim([0,max(histo)+2])
    xlabel('b [m]')
    ylabel('# solutions')
    title(['Average = ',num2str(vl(1)), ' m'])
    
    
    %%%%%%%%%%%%%%%%%%%
    W0H=(round(min(W0_final),3)-0.002  : 0.002  :  round(max(W0_final),3)+0.002);
    
    xbins=W0H(1:end-1);
    histo=[];
    for i=1:length(xbins)-1
        
        auxil=find(W0_final(:)<=xbins(i+1) & W0_final(:)>=xbins(i));
        rr=isempty(auxil) ;
        histo(i)=length(auxil);
    end
    
    bins=xbins(1:(end-1))+(0.0025/2);
    xi = linspace(min(bins), max(bins), 500);         
    yi = interp1(bins, histo, xi, 'v5cubic');
    subplot(4,2,6);
    hold off
    a=area(xi, yi);
    a.FaceAlpha = 0.25;
    %Hold on
    %plot(bins, histo)
    
    yl1=[0 30];
    vl=[mean(W0_final) mean(W0_final)];
    hold on
    plot(vl,yl1,'r--')
    
    xlim([min(bins)-0.002, max(bins)+0.002])
    ylim([0,max(histo)+2])
    xlabel('W_0 [m]')
    ylabel('# solutions')
    title(['Average = ',num2str(vl(1)), ' m/s'])
    
    a=findobj(gcf); % get the handles associated with the current figure
    allaxes=findall(a,'Type','axes');
    alllines=findall(a,'Type','line');
    alltext=findall(a,'Type','text');
%     set(allaxes,'FontName','Arial','LineWidth',1.2,...
%         'FontSize',11);
%     set(alllines,'Linewidth',1);
%     set(alltext,'FontName','Arial','FontSize',16);
    x0p=10;  y0p=220;
    widthp=1000;  heightp=1000;
    set(gcf,'position',[x0p,y0p,widthp,heightp])
    
    
    toc
    
    
    %%
    save('Results_Coppers','W0_final','B_final', 'T_final', 'L_final')
    
    % Distance I know
    X0=mean(X0_p);
    Y0=mean(Y0_p);
    if alpha<pi/2
        theta = alpha-atan(X0/Y0);
    else
        theta = atan(X0/Y0)-pi+alpha;
    end
    
    % Calculated by me
    X0_new;
    Y0_new;
    
    clc
    X = ['X0 calculated by the analytical solution: ',num2str(X0_new),' m.'];
    disp(X)
    X = ['Y0 calculated by the analytical solution: ',num2str(Y0_new),' m.'];
    disp(X)
    X = ['X0 retrieved by the model: ',num2str(X0),' m.'];
    disp(X)
    X = ['Y0 retrieved by the model: ',num2str(Y0),' m.'];
    disp(X)
    
    dd=sqrt(X0_new^2+Y0_new^2);
    X = ['Distance calculated by the analytical solution: ',num2str(dd),' m.'];
    disp(X)
    X = ['Real distance: ',num2str(r_dist),' m.'];
    disp(X)
    
    
    else
        
     
        hhh1=figure(1);
   
    hb1=subplot(4,2,5);
    hold off
    %title(['no tsunami potential'],'color','green')
    box on
    axis on
    
    xlabel('b [m]')
    ylabel('# solutions')
    
    
    hb2=subplot(4,2,6);
    
    %title(['no tsunami potential'],'color','green')
    box on
    axis on
    
    xlabel('W_0 [m]')
    ylabel('# solutions')
    
    hb3=subplot(4,2,7);
    hold off
    %title(['no tsunami potential'],'color','green')
    box on
    axis on
    xlabel('L [m]')
    ylabel('# solutions')
    hb4=subplot(4,2,8);
    hold off
    %title(['no tsunami potential'],'color','green')
    box on
    axis on
    xlabel('T [s]')
    ylabel('# solutions')
    
    
    a=findobj(gcf); % get the handles associated with the current figure
    allaxes=findall(a,'Type','axes');
    alllines=findall(a,'Type','line');
    alltext=findall(a,'Type','text');
%     set(allaxes,'FontName','Arial','LineWidth',1.2,...
%         'FontSize',11);
%     set(alllines,'Linewidth',1);
%     set(alltext,'FontName','Arial','FontSize',16);
    x0p=10;  y0p=220;
    widthp=1000;  heightp=1000;
    set(gcf,'position',[x0p,y0p,widthp,heightp])

    

end

% 
% figure(1)
% 
% frame = getframe(hhh1); 
%       im = frame2im(frame); 
%       [imind,cm] = rgb2ind(im,256); 
%       
%      if ii == 2 
%           imwrite(imind,cm,'signals.gif','gif','DelayTime',1,'Loopcount',inf); 
%       else 
%           imwrite(imind,cm,'signals.gif','gif','DelayTime',1,'WriteMode','append'); 
%      end 
%   

    
%end
b_m = mean(b_final);
L_m = mean(L_final);
T_m = mean(T_final);
W0_m = mean(W0_final);


