
function Direct_Model_2022_01_03(b, L, tau, W_0, X, Y, H,HOTSPOT,index)

for j=1:length(index)
%%
rhol=1000;%kg/m^3
Cl=1500;%m/s
%h=4000;
g=10;
%b=40000;%length width
%tau=10;%s
eta0=tau*W_0;%m
%x=1000000;
%y=0;
z=0;
n=10;
deltat=0.05;
%L=400000;%fault length

%%
%t0=(sqrt(X(j)^2+Y(j)^2))/sqrt(H(j)*g); 
t0=X(j)/sqrt(H(j)*g); 
t=t0:10:t0+2000;

%t0_a = (sqrt(X(j)^2+Y(j)^2))/Cl;
t0_a = X(j)/Cl;
t_a= t0_a:1:t0+1000;

%    for i=1:length(t) 
%     [ETA(i,1),PT(i,1),ETA0(i,1),Pb0(i,1)]=DM_stiassnie(H(j),b,tau,eta0,X(j),t(i),n);
%    end



[ETA0slender,ETAslender,pressure0,pressure]= DM_williams(H(j),b,L,tau,eta0,X(j),Y(j),H(j),t,n, t_a);                                                                       
% t0=(sqrt(X(j)^2+Y(j)^2))/sqrt(H(j)*g);
%t0=X(j)/sqrt(H(j)*g);
%[i0,j0]=find(abs(t-t0)==min(abs(t-t0)));
%ETA0slender(1:i0,1)=nan;
%%
width=1000;  % Width of figure for movie [pixels]
height=1000;  % Height of figure of movie [pixels]
left=20;     % Left margin between figure and screen edge [pixels]
bottom=20;  % Bottom margin between figure and screen edge [pixels]

%% plot surface elevation at HotSpots
figure(20)
set(gcf,'Position', [left bottom width height])

s2=subplot(length(index),1,j);

plot(t,real(ETA0slender));
hold on

axis on
box on
grid on
 x0p=1000;  y0p=600;
    widthp=750;  heightp=470;
    set(gcf,'position',[x0p,y0p,widthp,heightp])

%hold off

% if j==1
% l=legend('Williams et al. (2021)','Stiassnie (2010)','Orientation','horizontal','Location','southeast');
% l.EdgeColor='w';
% l.Box='off';
% l.FontSize=10;
% end

 if j==length(1)
          title(['$\eta$ [m]'],'interpreter','latex','Fontsize',16,'FontWeight','bold')
 end
 if j==length(index)
          xlabel('t [s]','interpreter','latex','Fontsize',10)
 end
 %if j~=length(index)
 %         set(gca,'xticklabel',[])
 %end

        ylabel(['#',HOTSPOT.hotspot_name{index(j)}],'interpreter','latex','Fontsize',10);

%       box on
%       axis on
%  
%       set(gca,'FontSize',10); 
      
      
      
%% plot pressure at HotSpots
figure(30)
set(gcf,'Position', [left 0 width height])
s2=subplot(length(index),1,j);
plot(t_a,real(pressure),'color','r');
hold on
%plot(t,PT,'color','k');

%ylim([-1000 1000])
%xlim([670 770])
axis on
box on
grid on
hold off
 x0p=1000;  y0p=0;
    widthp=750;  heightp=470;
    set(gcf,'position',[x0p,y0p,widthp,heightp])

% if j==1
% l=legend('Williams et al. (2021)','Stiassnie (2010)','Orientation','horizontal','Location','southeast');
% l.EdgeColor='w';
% l.Box='off';
% l.FontSize=10;
% end

 if j==length(1)
          title(['Pressure [Pa]'],'interpreter','latex','Fontsize',16,'FontWeight','bold')
 end
 if j==length(index)
          xlabel('t [s]','interpreter','latex','Fontsize',10)
 end
% if j~=length(index)
%          set(gca,'xticklabel',[])
% end

        ylabel(['#',HOTSPOT.hotspot_name{index(j)}],'interpreter','latex','Fontsize',10);

%       box on
%       axis on
%  
      set(gca,'FontSize',10); 
      
      
      
      
end

   %print(gcf,'-dpdf','Surface at 1000 km.pdf');