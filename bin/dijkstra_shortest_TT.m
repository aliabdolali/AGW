function [dist,path,x_transect,y_transect,z_transect,r_transect] = dijkstra_shortest_TT(x,y,depth,tri,Epi)

%Ali Abdolali EMC/NCEP/NOAA ali.abdolali@noaa.gov 22, March 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dijkstra_shortest_TT program Calculates the path that takes the shortest time between
% the epicenter and the rest of poist on the input mesh (triangular unstrucured)
% Some portion of the codes are taken from the original script written by
% Joseph Kirk but it is adopted for tsunami travel time, with shortest time
% via the fastest path in ocean using triangular unstrucured meshes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%inputs:
% x:      Longitude of nodes (deg)                           (N,1) 
% y:      Latitude of nodes (deg)                            (N,1) 
% depth:  Depth of nodes (m)                                 (N,1) 
% tri:    Element connections                                (M,3)
% Epi:    Longitude and latitude of the earthquake epicenter (1,2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs:
%     dist: The shortest travel time  (INF for any nodes that cannot be
%     reached along connection on the map).
%     path: ID list for the shortest path from epicenter to all the
%     nodes on the mesh (NAN for no direct path).
%     x_transect: the longiutde of the shortest path from epiceter to all
%     nodes (deg)
%     y_transect: the latiutde of the shortest path from epiceter to all
%     nodes (deg)
%     z_transect: the depth of the shortest path from epiceter to all nodes
%     (m)
%     r_transect: the cumulative  distance of the shortest path from
%     epiceter to all nodes (km)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nodes(:,1)=1:length(x);
nodes(:,2)=x;
nodes(:,3)=y;

%Put together all the connection, sort and remove the repeated ones. 
m=length(tri(:,1));
s(:,1)=tri(:,1);s(:,2)=tri(:,2);
s2(:,1)=tri(:,1);s2(:,2)=tri(:,3);
s3(:,1)=tri(:,2);s3(:,2)=tri(:,3);
s(m+1:2*m,1:2)=s2;s(2*m+1:3*m,1:2)=s3;
ss=sortrows(s);
sss = unique(ss(:,[1 2]),'rows');
segments(:,1)=1:length(sss(:,1));
segments(:,2:3)=sss;

%find the nearest node ID to the epicenter
[start_id,DIS]=knnsearch([x,y],[Epi(1),Epi(2)],'k',1);
% 

    % initializations
    node_ids = nodes(:,1);
    [num_map_pts,cols] = size(nodes);
    table = sparse(num_map_pts,2);
    shortest_distance = Inf(num_map_pts,1);
    settled = zeros(num_map_pts,1);
    path = num2cell(NaN(num_map_pts,1));
    col = 2;
    pidx = find(start_id == node_ids);
    shortest_distance(pidx) = 0;
    table(pidx,col) = 0;
    settled(pidx) = 1;
    path(pidx) = {start_id};
 
        while_cmd = 'sum(~settled) > 0';

    while eval(while_cmd)
        % update the table
        table(:,col-1) = table(:,col);
        table(pidx,col) = 0;
        % find neighboring nodes in the segments list
        neighbor_ids = [segments(node_ids(pidx) == segments(:,2),3);
            segments(node_ids(pidx) == segments(:,3),2)];
        % calculate the distances to the neighboring nodes and keep track of the paths
        for k = 1:length(neighbor_ids)
            cidx = find(neighbor_ids(k) == node_ids);
            if ~settled(cidx)
                %d = sqrt(sum((nodes(pidx,2:cols) - nodes(cidx,2:cols)).^2));
                %[d1km d2km]=lldistkm([nodes(pidx,3) nodes(pidx,2)],[nodes(cidx,3) nodes(cidx,2)]);
                [d1]=latlondistm([nodes(pidx,3) nodes(pidx,2)],[nodes(cidx,3) nodes(cidx,2)]);
                %load dA
                %depth=dA;
                d=d1./sqrt(9.81*depth(pidx));
                %d=d1;
                %d=1000*d1km;
                if (table(cidx,col-1) == 0) || ...
                        (table(cidx,col-1) > (table(pidx,col-1) + d))
                    table(cidx,col) = table(pidx,col-1) + d;
                    tmp_path = path(pidx);
                    path(cidx) = {[tmp_path{1} neighbor_ids(k)]};
                else
                    table(cidx,col) = table(cidx,col-1);
                end
            end
        end
        % find the minimum non-zero value in the table and save it
        nidx = find(table(:,col));
        ndx = find(table(nidx,col) == min(table(nidx,col)));
        if isempty(ndx)
            break
        else
            pidx = nidx(ndx(1));
            shortest_distance(pidx) = table(pidx,col);
            settled(pidx) = 1;
        end
    end

        dist = shortest_distance';
        path = path';
        % extract the depth
        for i=1:length(path)
            if ~isnan(sum(path{i}))
               x_transect{i}=x(path{i});
               y_transect{i}=y(path{i});
               z_transect{i}=depth(path{i});
               r_transect{i}=[0; cumsum(latlondistm([y(path{i}(2:end)) x(path{i}(2:end))]...
                                        ,[y(path{i}(1:end-1)) x(path{i}(1:end-1))]))];
            else
            x_transect{i}=nan; 
            y_transect{i}=nan;  
            z_transect{i}=nan; 
            r_transect{i}=nan; 
            
            end
        end
       
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