clear all
clc
%This script reads ETOPO1, and draw an arc between the epicenter and a dart
%buoy.
%Ali Abdolali EMC/NCEP/NOAA ali.abdolali@noaa.gov 22, March 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input
Epi=[142.369,38.322]; %coordinare of the epicenter
dx=0.1;       % Transect resolution
I=24; %Index of the NDBC
new_origin=[0,-220]; %for changing the origin of plot (the default is [0 0] 
%which plots [-180 180]; new_origin=[0,-220] is for source and destiation
%in the pacific
debug=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load etopo topobathy
ncfile='etopo1.nc';
lon1=ncread(ncfile,'lon')-180;
[ip,jp]=find(lon1>=0);
[in,jn]=find(lon1<0);
lat=ncread(ncfile,'lat');
depth1=flipud(rot90(ncread(ncfile,'z')));
depth(:,1:length(ip))=depth1(:,ip);
depth(:,end+1:end+length(in))=depth1(:,in);
lon(1:length(ip),1)=lon1(ip)-180;
lon(end+1:end+length(in),1)=lon1(in)+180+(lon1(2)-lon1(1));
[LON,LAT]=meshgrid(lon,lat);
clear lon1
clear depth1
display('Read ETOPO1 ...');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read the DART coordinates
S=shaperead('NDBC.shp');
for i=1:length(S)
NDBC{i,1}=S(i).name;
xNDBC(i,1)=S(i).X;
yNDBC(i,1)=S(i).Y;
end

display(['EPICENTER = [',num2str(Epi(1)),' ',num2str(Epi(2)),']']);
display(['NDBC ',NDBC{I}]);
display(['NDBC Coordinates = [',num2str(xNDBC(I)),' ',num2str(yNDBC(I)),']']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[arclen,az] = distance('gc',[Epi(2),Epi(1)],[yNDBC(I),xNDBC(I)]); %[arc length, azimuth]
npnt=floor(arclen/dx); %number of points on the arc
waypoints=[Epi(2),Epi(1);yNDBC(i),xNDBC(i)];
[lattrkgc,lontrkgc] = track2('gc',Epi(2),Epi(1),yNDBC(I),xNDBC(I),[1 0],'degrees',npnt);
%az = azimuth('gc',Epi(2),Epi(1),yNDBC(I),xNDBC(I))
[lttrk,lntrk] = track('gc',waypoints,'degrees');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract depth along arc
Rtransect(1)=0;
Ztransect=interp2(LON,LAT,depth,lontrkgc,lattrkgc);
for i=1:length(lontrkgc)-1
    [d1km(i+1) d2km(i+1)]=lldistkm([lattrkgc(i+1),lontrkgc(i+1)],[lattrkgc(i),lontrkgc(i)]);
Rtransect(i+1)=Rtransect(i)+d1km(i+1);
end
display('EXTRACT TRANSECT ...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load coast
if debug==1
figure
pcolor(LON(1:100:end,1:100:end),LAT(1:100:end,1:100:end),depth(1:100:end,1:100:end));
shading interp
colormap('jet')
view(2)
hold on
scatter(xNDBC,yNDBC,'ok','filled')
hold on
scatter(Epi(1),Epi(2),'sb','filled')
hold on
scatter(xNDBC(I),yNDBC(I),'sb','filled')
hold on
scatter(lontrkgc,lattrkgc,'xb')
hold on
scatter(lntrk,lttrk,'xr')


figure
load coastlines
axesm eckert4
setm(gca,'mapprojection','ortho')
geoshow(LAT(1:100:end,1:100:end),LON(1:100:end,1:100:end),depth(1:100:end,1:100:end),'DisplayType','surface')
framem
gridm
plotm(Epi(2),Epi(1),'rp')
hold on
plotm(yNDBC,xNDBC,'m.')
%axesm('mercator','MapLatLimit',[-90 90],'MapLonLimit',[-180 180],...
%    'Grid','on','Frame','on','MeridianLabel','on','ParallelLabel','on')
geoshow(coastlat,coastlon,'DisplayType','line','color','b')
hold on
geoshow(lattrkgc,lontrkgc,'DisplayType','line','color','r')
%framem
%gridm

end
%%

load coastlines

width=680;  % Width of figure for movie [pixels]
height=1000;  % Height of figure of movie [pixels]
left=200;     % Left margin between figure and screen edge [pixels]
bottom=200;  % Bottom margin between figure and screen edge [pixels]

  figure
set(gcf,'Position', [left bottom width height])
s1= subplot(2,1,1)  
[EpiR(2),EpiR(1)]=rotatem(Epi(2),Epi(1),new_origin,'inverse','degrees');
[yNDBCR,xNDBCR]=rotatem(yNDBC,xNDBC,new_origin,'inverse','degrees');
[coastlatR,coastlonR]=rotatem(coastlat,coastlon,new_origin,'inverse','degrees');
[lattrkgcR,lontrkgcR]=rotatem(lattrkgc,lontrkgc,new_origin,'inverse','degrees');

axesm eckert4
%%% if you want the plot on a sphere
%setm(gca,'mapprojection','ortho') 
%geoshow(LAT(1:100:end,1:100:end),LON(1:100:end,1:100:end),depth(1:100:end,1:100:end),'DisplayType','surface')
framem
gridm
plotm(EpiR(2),EpiR(1),'rp')
hold on
plotm(yNDBCR,xNDBCR,'m.')
%axesm('mercator','MapLatLimit',[-90 90],'MapLonLimit',[-180 180],...
%    'Grid','on','Frame','on','MeridianLabel','on','ParallelLabel','on')
geoshow(coastlatR,coastlonR,'DisplayType','line','color','b')
hold on
geoshow(lattrkgcR,lontrkgcR,'DisplayType','line','color','r')
hold on
%framem
%gridm

title(['EPICENTER [',num2str(Epi(1)),' ',num2str(Epi(2)),'] - ','DART#',NDBC{I}])
     
s2=subplot(2,1,2)
clear poly
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
ylim([min(Ztransect)/1000-1,1])
xlim([-100,Rtransect(end)+100])


e=scatter(Rtransect(1),Ztransect(1)/1000,80,'p')
e.MarkerEdgeColor = 'k';
e.MarkerFaceColor = 'r';

s=scatter(Rtransect(end),0,60,'v')
s.MarkerEdgeColor = 'm';
s.MarkerFaceColor = 'm';


l=legend([e,s],'Epicenter',['DART#',NDBC{I}],'Orientation','horizontal','Location','southeast')
l.EdgeColor='w'
l.Box='off'
l.FontSize=12

         xlabel('Distance [km]','interpreter','latex','Fontsize', 8)
        ylabel('h [km]','interpreter','latex','Fontsize', 8);
           
    
      box on
      axis on
 
      set(gca,'FontSize',12); 
      
      
      
      set(s1, 'Position', [.08 .52 .88 .4]);
      set(s2, 'Position', [.08 .06 .88 .4]);
      
print(gcf,'-dpng',['DART',NDBC{I},'.png'],'-r500');
