clear all
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program downloads the hidtorical DART data for a given year       %
% Ali Abdolali (EMC/NCEP/NOAA ali.abdolali@noaa.gov                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%    INPUT    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yr=2022;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B={'21346' '21347' '21348' '21401' '21402' '21413' '21414' '21415' '21416' '21417' '21418' '21419'...
'21420' '23217' '23218' '23219' '23220' '23223' '23226' '23227' '23228' '23401' '23461' '32066'...
'32067' '32068' '32069' '32401' '32402' '32403' '32404' '32411' '32412' '32413' '32489' '34420'...
'41420' '41421' '41424' '41425' '42407' '42408' '42409' '43412' '43413' '44401' '44402' '44403'...
'46010' '46401' '46402' '46403' '46404' '46405' '46406' '46407' '46408' '46409' '46410' '46411'...
'46412' '46413' '46414' '46415' '46416' '46419' '46451' '46452' '46g10' '51406' '51407' '51425'...
'51426' '52401' '52402' '52403' '52404' '52405' '52406' '53046' '53401' '54401' '55012' '55013'...
'55015' '55016' '55023' '55042' '55401' '56001' '56003' 'dartp' 'dartq' 'dartr' 'darts' 'dartt'...
'dartu'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%swden - Spectral Wave Density data with Spectral Wave Direction data
disp('DART Buoy Data')
disp('Start ...')
disp(['year = ',num2str(yr)])
if num2str(yr)==datestr(date,'yyyy')
    disp(['close to real time'])
end
m=0; 


for i=1:length(B)
    %check if it exist locally
    if isfile(['../reference_data/',B{i},'t',num2str(yr),'.nc'])
       disp([B{i},'t',num2str(yr),'.nc exists.'])
    else
myURL=['https://dods.ndbc.noaa.gov/thredds/fileServer/data/dart/',B{i},'/',B{i},'t',num2str(yr),'.nc'];
myURL9999=['https://dods.ndbc.noaa.gov/thredds/fileServer/data/dart/',B{i},'/',B{i},'t9999.nc'];

[str,status] = urlread(myURL);
%if wave data is available
if status==1
    m=m+1;
urlwrite(myURL,['../reference_data',B{i},'t',num2str(yr),'.nc']);
%XY1(m,1)=X(i);
%XY1(m,2)=Y(i);
II1(m,1)=i;
disp(['Downloading DART#',B{i},' ...']) 
elseif num2str(yr)==datestr(date,'yyyy')
    [str,status] = urlread(myURL9999);
    %if wave data is available
      if status==1
        m=m+1;
      urlwrite(myURL9999,['../reference_data/',B{i},'t',num2str(yr),'.nc']);
      %XY1(m,1)=X(i);
      %XY1(m,2)=Y(i);
      II1(m,1)=i;
       disp(['Downloading DART#',B{i},' ...'])
      end
else
disp(['DART#',B{i},' is not avaiable']) 
end
end
end
disp('Finished.')