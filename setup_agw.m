function setup_agw
% -------------------------------------------------------------------------
%|                                                                        |
%|                    +----------------------------+                      |
%|                    | setup_agw                  |                      |
%|                    |                            |                      |
%|                    | Last Update: 21-March-2021 |                      |
%|                    +----------------------------+                      |
%|                     Distributed with AGW package                       |
%|                                                                        |
%|                 Copyright ................,                            |
%|                     All rights reserved.                               |
%| Ali Abdolali (ali.abdolali@noaa.gov & Usama Kadri (kadriu@cardiff.ac.uk| 
%|                                                                        |
%| DESCRIPTION                                                            |
%| The setup_agw function supports the Acoustic Gravity Waves Exapmles    |
%| 1. downloads the reference data (etopo1 and etopo2) from NCEP's server,|
%| 2. uncompresses the tarball at the appropriate path and                |
%| 3. adds temporarily the agw_tools's paths to the MATLAB's pathdef.     |
%|                                                                        |
%| INPUT                                                                  |
%| None                                                                   | 
%|                                                                        |
%| OUTPUT                                                                 |
%| None                                                                   |
%|                                                                        |
%| NOTES                                                                  |
%| In case of updates update the section "Define paths, files and ftp     |
%| server and paths".                                                     |
%|                                                                        |
%| BUG FIXES                                                              |
%|                                                                        |
%| PRGRMR   :   Ali Abdolali                                              |
%| DATE     :                                                             |
%|              v.1.0 - 21-March-2021                                     |
% -------------------------------------------------------------------------
%%
display('AGW_tools installation!')
%% Define paths, files and ftp server and paths
home = fileparts(which(mfilename)); % AGW_tools directory
path_tar='reference_data';          % reference data path
path_bin='bin';                     % bin path
path_exm='examples';                % examples path
%
ftp_svr='polar.ncep.noaa.gov';      % ftp server of reference data
ftp_pth='/tempor/ww3ftp';           % ftp path for reference data
bathy_file='gridgen_addit.tgz';     % reference data tarball
%% Downloading
if exist([home,'/',path_tar,'/etopo1.nc'], 'file') ~= 2
    ftp_ind=ftp(ftp_svr);
    cd(ftp_ind,ftp_pth);
    mget(ftp_ind,bathy_file,home);
    close(ftp_ind);
%% Untar the reference data
    untar([home,'/',bathy_file],path_tar);
%
    delete([home,'/',bathy_file]);
end
%download unstr mesh for shortest path calculation
if exist([home,'/',path_tar,'/global60_50km_unstr.msh'], 'file') ~= 2
    urlwrite('https://drive.google.com/file/d/1OKC6aHnTf-lKwncOwBJvQ2jsxsuv2MBt/view?usp=sharing',[home,'/',path_tar,'/global60_50km_unstr.msh']);
end
%% Add the bin, reference path and examples to the user's matlab path
addpath(fullfile(home, path_bin))
addpath(fullfile(home, path_tar));
addpath(fullfile(home, path_exm));

%% Install NetCDF package for octave
vers=ver;
% interpreter='Matlab';
% netcdf_install='false';
for i1=1:1:length(vers)
    if strcmpi (vers(i1).Name, 'Octave')
        pkg install -forge netcdf
    end
end
