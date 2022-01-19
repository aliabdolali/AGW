function [OUT]=hotspot_model_20220107(Epi,hotspot,ttt)
%%
%clear all
%clc

%input
%epicenter coordinate
%Epi=[143.76,27.10];
%name of hotspots
%N={'   21413' '   51425' '   51426' '   54401' '   21346' '   21347' '   21401'};
%latitude of hotspots
%yN=[30.517; -9.505;    -23.110;   -33.109; 40.302;  39.601;   42.617];
%longitude of hotspots
%xN=[152.127; -176.262; -168.385; -173.155; 146.192; 145.799; 152.583];
%plot tsunami travel time if ttt=true
%%
%output
%OUT.depth_Epi; epicenter depth
%OUT.Zmean; mean of each transect from epicenter to hotspot
%OUT.az; azimuth between epicenter and hotspot
%OUT.arclendeg; Arc length from epicenter to hotspot (degree)
%OUT.arclenkm; Arc length from epicenter to hotspot (km)
%OUT.transect_latitude; shortest path latitude
%OUT.transect_longitude; shortest path longitude
%OUT.transect_distance_km ; shortest path distance
%OUT.transect_depth_m; shortest path depth profile
%%

xN=hotspot.xN;
yN=hotspot.yN;
N=hotspot.N;
%%
%read depth file (trinagular mesh)_
I=1:length(xN);
npnt=500;
[triA,xA,yA,dA] = readWW3mesh('global60_50km_unstr.msh',0);
dA(dA<=10)=0;
%%
%calculate shortest path
[DIST,path,x_transect,y_transect,z_transect,r_transect] = dijkstra_shortest_TT(xA,yA,dA,triA,Epi);
%calculate epicenter depth
depth_Epi=griddata(xA,yA,dA,Epi(1),Epi(2));
%phase speed of gravity wave (shallow water approximation)
Gravity=sqrt(9.81*dA);

if strcmp(ttt,'true')
%interpolate on a regular grid (curvilinear)
xR=-180:0.5:180;
yR=-85:0.5:85;
[XR,YR]=meshgrid(xR,yR);
DISTR=griddata(xA,yA, DIST,XR,YR);
end
%%
if strcmp(ttt,'true')
%calculate contours for the plot
clear poly
clear polx
clear labell
m=1;
mm=0;
[c1,h1]=contour(YR,XR,DISTR,[0:300:15*3600]);
while m<=length(c1(1,:)) 
  poly{mm+1}=c1(1,m+1:m+c1(2,m));  
  polx{mm+1}=c1(2,m+1:m+c1(2,m));  
  labell(mm+1)=c1(1,m);
  mm=mm+1;
  m=m+c1(2,m)+1;
end
end
close all
%%
%calculate the shortest arc
dest=[xN(I) yN(I)];
for i=1:length(I)
[lattrkgc(:,i),lontrkgc(:,i)] = track2('gc',Epi(2),Epi(1),dest(i,2),dest(i,1),[1 0],'degrees',npnt);
[arclendeg(i,1),az(i,1)] = distance('gc',Epi(2),Epi(1),dest(i,2),dest(i,1));
arclenkm(i,1)=deg2km(arclendeg(i,1));
%[latT(i),lonT(i)] = track1(Epi(2),Epi(1),az(i),arclen(i)) 
%X = cosd(dest(i,2))* sind(dest(i,1)-Epi(1));
%Y = cosd(Epi(2)) * sind(dest(i,2)) - sind(Epi(2)) * cos(dest(i,2)) * cosd(dest(i,1)-Epi(1));
%beta(i) = atan2(X,Y)*180/pi;
end

%r=zeros(size(I));
%rr=r;
%for i=1:length(I)
%  for j=1:npnt-1
%   [d1]=latlondistm([lattrkgc(j,i) lontrkgc(j,i)],[lattrkgc(j+1,i) lontrkgc(j+1,i)]);
%    r(j,i)=d1;
%    rr(i)=rr(i)+d1;
%  end
%end
%%
%find the node in the DIST for shortest path 
[sid,d]=knnsearch([xA,yA],[Epi(1),Epi(2)],'k',1);
[fid,d]=knnsearch([xA,yA],[dest(:,1),dest(:,2)],'k',1);

for i=1:length(I)
Zmean(i)=nanmean(z_transect{fid(i)});
end

%%
%output
OUT.depth_Epi=depth_Epi;
OUT.Epi=Epi;
OUT.Zmean=Zmean;
OUT.az=az;
OUT.arclendeg=arclendeg;
OUT.arclenkm=arclenkm;
OUT.hotspot_longitude=hotspot.xN;
OUT.hotspot_latitude=hotspot.yN;
OUT.hotspot_name=hotspot.N;
for i=1:length(I)
OUT.transect_latitude{i}=yA(path{fid(i)});
OUT.transect_longitude{i}=xA(path{fid(i)});
OUT.transect_distance_km{i}=r_transect{fid(i)};
OUT.transect_depth_m{i}=z_transect{fid(i)};
end
%%
if strcmp(ttt,'true')
DISTRHR=DISTR/3600;
DISTRHR(:,1)=DISTRHR(:,2);
DISTRHR(:,end)=DISTRHR(:,end-1);
hmax=15;
clear mymap
mymap(:,1)=0:1/((2*hmax)-1):1;
mymap(:,2)=0:1/((2*hmax)-1):1;
mymap(:,3)=0:1/((2*hmax)-1):1;
end
width=1000;  % Width of figure for movie [pixels]
height=1000;  % Height of figure of movie [pixels]
left=20;     % Left margin between figure and screen edge [pixels]
bottom=20;  % Bottom margin between figure and screen edge [pixels]
figure
set(gcf,'Position', [left bottom width height])

%latlim = [nanmin(Epi(2),nanmin(yN))-10 nanmin(Epi(2),nanmin(yN))+10];
%lonlim = [nanmin(Epi(1),nanmin(xN))-10 nanmin(Epi(1),nanmin(xN))+10];

latlim = [-50 50];
lonlim = [100 250];


axesm eckert4
%axesm ortho
axesm('robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,...
    'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on','MLabelParallel','south','FontSize',18)
if strcmp(ttt,'true')
caxis([0 hmax])%this one is user defined
levels=[0:1:24];
hold on
geoshow(YR,XR,DISTRHR,'DisplayType','texturemap')
hold on
end
geoshow(yA(fid),xA(fid),'DisplayType','multipoint','Marker', 'o','MarkerSize',8, 'Color', 'r', 'MarkerFaceColor', 'r','MarkerEdgeColor', 'k')
hold on
geoshow(yA(sid),xA(sid),'DisplayType','point', 'Marker', 'h','MarkerSize',20, 'Color', 'k', 'MarkerFaceColor', 'y','MarkerEdgeColor', 'auto')
hold on
textm(yA(fid),xA(fid)+0.3, N, 'FontSize',8,'color','b')
hold on

for i=1:length(I)
geoshow(yA(path{fid(i)}),xA(path{fid(i)}),'DisplayType','line','Color', 'm')
hold on
geoshow(lattrkgc(:,i),lontrkgc(:,i),'DisplayType','line','Color', 'r')
hold on
end
geoshow('landareas.shp','FaceColor',[0.2 0 0],'FaceAlpha', 1);
hold on
setm(gca,'MLabelLocation',60)
hold on
if strcmp(ttt,'true')
colormap(mymap)
hcb=colorbar

hcb.FontSize=18
hcb.LineWidth = 1.5;
l=[0:1:hmax];
ll=[0:1:hmax];
set(hcb,'XTick',l,'XTickLabel', ll);
title(hcb,'[hr]','Fontsize',18)
title(['Tsunami Travel Time (TTT)'],'fontsize',18)
end
box off
axis off
print(gcf,'-dpng',['TT_gravity_n','.png'],'-r900');


%%

width=1000;  % Width of figure for movie [pixels]
height=1000;  % Height of figure of movie [pixels]
left=20;     % Left margin between figure and screen edge [pixels]
bottom=20;  % Bottom margin between figure and screen edge [pixels]
figure
set(gcf,'Position', [left bottom width height])
for i=1:length(I)
R(i)=max(r_transect{fid(i)});
Z(i)=min(-z_transect{fid(i)});
end


for i=1:length(I)

s2=subplot(length(I),1,i)
clear poly
Rtransect=r_transect{fid(i)};
Ztransect=-z_transect{fid(i)};
poly(:,1)=Rtransect;
poly(:,2)=Ztransect;
poly(length(Rtransect)+1,1)=Rtransect(end);
poly(length(Rtransect)+1,2)=0;
poly(length(Rtransect)+2,1)=0;
poly(length(Rtransect)+2,2)=0;


h1=fill(poly(:,1),poly(:,2)/1000,[135/255 206/255 250/255]);
set(h1,'edgecolor','white');

hold on
plot(Rtransect,Ztransect/1000,'k','linewidth',1)
hold on
plot(Rtransect,0*Ztransect/1000,'k','linewidth',1)
hold on
plot([0 0],[Ztransect(1)/1000 0],'--k','linewidth',2)
hold on
plot([Rtransect(end) Rtransect(end)],[Ztransect(end)/1000 0],'--k','linewidth',2)
hold on

xlim([-100,max(R(:))+100])
ylim([floor(min(Z)/1000)-1,1])

e=scatter(Rtransect(1),Ztransect(1)/1000,80,'p')
e.MarkerEdgeColor = 'k';
e.MarkerFaceColor = 'r';

s=scatter(Rtransect(end),0,60,'v')
s.MarkerEdgeColor = 'm';
s.MarkerFaceColor = 'm';

text(10,floor(min(Z)/1000),['#',N{i}],'Fontsize', 14)
if i==1
l=legend([e,s],'Epicenter','Obs','Orientation','horizontal','Location','southeast')
l.EdgeColor='w'
l.Box='off'
l.FontSize=18
end

if i==I(end)
         xlabel('Distance [km]','interpreter','latex','Fontsize', 18)
end
if i~=I(end)
         set(gca,'xticklabel',[])
end

        ylabel(['h [km]'],'interpreter','latex','Fontsize', 168);

      box on
      axis on
 
      set(gca,'FontSize',18); 
      
end
      
print(gcf,'-dpng',['transects_gravity_n','.png'],'-r900');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d1]=latlondistm(latlon1,latlon2)
radius=6371;
   lat1=latlon1(:,1)*pi/180;
   lat2=latlon2(:,1)*pi/180;
   lon1=latlon1(:,2)*pi/180;
   lon2=latlon2(:,2)*pi/180;
   deltaLat=lat2-lat1;
   deltaLon=lon2-lon1;
   a=sin((deltaLat)/2).^2 + cos(lat1).*cos(lat2).*sin(deltaLon/2).^2;
   c=2*atan2(sqrt(a),sqrt(1-a));
   d1(:,1)=1000*radius*c;    %Haversine distance
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%