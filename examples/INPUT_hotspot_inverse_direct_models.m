% Inverse problem model input

close all
clear 
clc 

%inputs

min_f=1;
%min_f= minimum possible 1st mode frequency
max_f=10;
%max_f= maximum possible 1st mode frequency
npts=4;
%npts= input(prompt) % NUMBER OF POINTS FROM THE SIGNAL THAT I WANT FOR COMPUTING EACH SET OF SOLUTIONS
nnn=10;
%nnn = input(prompt) % NUMBER OF SETS OF SOLUTIONS THAT I WANT
Epi=[143.76,27.10];
%epicenter coordinate (lon,lat)
hotspot.N={ '  51425'   '  21401'   '  HA11' '  Hotspot1' '  Hotspot2' '  Hotspot3'};
%name of hotspots
hotspot.xN=[ -176.262; 152.583; 166.6; 135; -160; -125];
%longitude of hotspots
hotspot.yN=[-9.505;  42.617; 19.3; 31; 19.3; 37];
%latitude of hotspots
ttt='false';
%if ttt=='true', the code generates tsunami travel time plot as well
format longEng
%% execute hotspot model to estimate the distance and average depth
[HOTSPOT]=hotspot_model(Epi,hotspot,ttt);
%HOTSPOT.depth_Epi; epicenter depth
%HOTSPOT.Zmean; mean of each transect from epicenter to hotspot
%HOTSPOT.az; azimuth between epicenter and hotspot
%HOTSPOT.arclendeg; Arc length from epicenter to hotspot (degree)
%HOTSPOT.arclenkm; Arc length from epicenter to hotspot (km)
%HOTSPOT.transect_latitude; shortest path latitude
%HOTSPOT.transect_longitude; shortest path longitude
%HOTSPOT.transect_distance_km ; shortest path distance
%HOTSPOT.transect_depth_m; shortest path depth profile
%%
%% NECESSARY INPUT VARIABLES
Fs=250; % Sampling frequency
I_Hydro=3;
%the index of hydrophone in the hotspot model
h=HOTSPOT.Zmean(I_Hydro); % Average water depth between hydrophone and epicentre
r_dist=HOTSPOT.arclenkm(I_Hydro)*10^3; % Distance between epicentre and hydrophone
alpha = HOTSPOT.az(I_Hydro)*pi/180;   % <<<<<<<<<<<<<<<<<<< INPUT from HS Model this is the angle between EQ and the HA calculated from the North clockwise
DT=100; % update the script
Pmax=30; %max value within time window

% Use if input manually class and magnitude
clas=1; 
magnitude=7.07;
%% READ THE FILE AND SET SPECIFICATIONS OF HYDROPHONE
current_f=pwd;
cd(['reference_data/MATLAB_signal/'])
x = dir('mag_74_21_12_2010.txt');
fileID = fopen(x.name,'r');

formatSpec = '%f';
p = fscanf(fileID,formatSpec);
cd(current_f)
aux=1:length(p); TT=aux./Fs; % Generate time vector related to the signal
p=p-mean(p); % Remove mean
p=p./(10^6); % Conversion from Micropascals to Pascals! %%%%%%%%%%%%%%



%% call inverse problem model function: 
for i=2:floor(TT(end)/DT)
   [b, L, T, W_0, theta, X0, Y0] = Inverse_Model(min_f, max_f, npts, nnn, Fs, h, r_dist, DT, Pmax, clas, magnitude, TT, p, i, alpha);
    
    % the output of Inverse_Model should be fed into the direct problem
%     % the two should run in parallel 
   if isnan(b)==0 && isnan(L)==0 && isnan(T)==0 && isnan(W_0)==0
    I_hotspot=[1 2 3 4 5 6]; 
    for j=1:length(I_hotspot)    
       X(I_hotspot(j)) = abs(HOTSPOT.arclenkm(I_hotspot(j))*1000* sin(HOTSPOT.az(I_hotspot(j))*pi/180-theta));  
       Y(I_hotspot(j)) = HOTSPOT.arclenkm(I_hotspot(j))*1000* cos(HOTSPOT.az(I_hotspot(j))*pi/180-theta);
    end
     Direct_Model(b, L, T, W_0,X,Y,HOTSPOT.Zmean,HOTSPOT,I_hotspot);
   end
end

