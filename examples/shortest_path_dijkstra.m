clear all
clc
%This script reads a triangular mesh and use dijkstra technique to find the 
% fastest path from the the epicenter and a given DART
%buoy.
%Ali Abdolali EMC/NCEP/NOAA ali.abdolali@noaa.gov 22, March 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
% number of points for track2
npnt=200;
% read the coordinates of DART network
S=shaperead('NDBC.shp');
for i=1:length(S)
NDBC{i,1}=S(i).name;
xNDBC(i,1)=S(i).X;
yNDBC(i,1)=S(i).Y;
end
I=24; %Index of the NDBC

% read the triangular mesh (node, depth, connections).
[tri,x,y,d] = readWW3mesh('global60_50km_unstr.msh',0);
% remove the overland topo
depth=d;
d(d<=10)=nan;
depth(x>178.5)=nan;
depth(x<=-178.5)=nan;
% epicenter coordinate (Tohoku 2011)
Epi=[142.369,38.322];
%calculate the travel time along the fastest route
[DIST,path,x_transect,y_transect,z_transect,r_transect] = dijkstra_shortest_TT(x,y,d,tri,Epi);

% extraxct the depth from epicenter to nodes


%%
dest=[xNDBC(I) yNDBC(I)];
% find the shortest arc between epicenter and hot points (NDBCs) 
for i=1:length(I)
[lattrkgc(:,i),lontrkgc(:,i)] = track2('gc',Epi(2),Epi(1),dest(i,2),dest(i,1),[1 0],'degrees',npnt);
end
% find the nearest nodes in the triangular mesh to the epicenter.
[sid,distance1]=knnsearch([x,y],[Epi(1),Epi(2)],'k',1);
% find the nearest nodes in the triangular mesh to the hot points.
[fid,distance2]=knnsearch([x,y],[dest(:,1),dest(:,2)],'k',1);

toc

%%
dist=DIST;
dist(x>178.5)=nan;
dist(x<=-178.5)=nan;
DISTT(:,1)=dist(1,:);


width=680;  % Width of figure for movie [pixels]
height=1000;  % Height of figure of movie [pixels]
left=200;     % Left margin between figure and screen edge [pixels]
bottom=200;  % Bottom margin between figure and screen edge [pixels]

  figure
set(gcf,'Position', [left bottom width height])
s1= subplot(2,1,1)  
trisurf(tri,x,y,DISTT/3600);
caxis([0 24])
shading interp
view(2)
hcb=colorbar
colormap(jet)

hold on
plot(x(sid), y(sid),'ro');
hold on
plot(x(fid), y(fid),'bo');
hold on

for i=1:length(I)
plot3(x(path{fid(i)}), y(path{fid(i)}),100*ones(size(path{fid(i)})),'.k')
hold on
plot3(lontrkgc(:,i), lattrkgc(:,i),100*ones(length(lattrkgc(:,i)),1),'.m')
end
xlim([-180 180])
ylim([-90 90])
xlabel('Longitude [ deg. ]')
ylabel('Latitude [ deg. ]')

hcb.Label.String = 'Travel Time (hr)';
title(['EPICENTER [',num2str(Epi(1)),' ',num2str(Epi(2)),'] - ','DART#',NDBC{I}])
     
s2=subplot(2,1,2)
clear poly
Rtransect=r_transect{fid};
Ztransect=-z_transect{fid};
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
      
      
      
      set(s1, 'Position', [.08 .52 .80 .4]);
      set(s2, 'Position', [.08 .06 .88 .4]);
      
print(gcf,'-dpng',['DART',NDBC{I},'.png'],'-r500');

