function [output] = transform_filter(dt,y,varargin)
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Low pass
% varargin{1} = frequency of "input";
% varargin{2} = type of filter (insert string'pb' for lowpass);
%
% High Pass
% varargin{1} = frequency of "input";
% varargin{2} = type of filter (insert string'pa' for highpass);
% Bandpass
% varargin{1} = 1st frequency of "input";
% varargin{2} = last frequency;
% varargin{3} = type of filter (insert string'bandpass' for bandpass);
% 
% output.ytfilt = trasformation of  the Fourier filter;
% output.yant = filtered signal;
% 
% -------------------------------------------------------------------------
%                           -> Ali Abdolali

n=length(y);
output.yt=fft(y);
fappo=0 : 1/n/dt : 1/dt;
output.fre(1:ceil(n/2))=fappo(1:ceil(n/2));
ifr=0;
for indice=n:-1:ceil(n/2)+1
    ifr=ifr+1;
    fchemis=fappo(ifr+1);
    output.fre(indice)=fchemis;
end

if length(varargin) == 2
    filtro = varargin{2};
    f_cut_off = varargin{1};
    
    if     filtro == 'pb'
        disp('------------------');
        disp('Low pass filter');
        disp('------------------');
        indici_dirt = find(output.fre > f_cut_off);
        output.ytfilt = output.yt;
        output.ytfilt(indici_dirt) = 0;
        output.yant = ifft(output.ytfilt);
    elseif filtro == 'pa'
        disp('------------------');
        disp('High pass filter');
        disp('------------------');
        indici_dirt = find(output.fre < f_cut_off);
        output.ytfilt = output.yt;
        output.ytfilt(indici_dirt) = 0;
        output.yant = ifft(output.ytfilt);
    end
    
elseif length(varargin) == 3
    filtro = varargin{3};
    f_cut_in = varargin{1};    
    f_cut_off = varargin{2};        
        if   filtro == 'bandpass'
        disp('------------------');
        disp('Band pass filter');
        disp('------------------');
        indici_dirt = find(output.fre < f_cut_in | output.fre > f_cut_off);
        output.ytfilt = output.yt;
        output.ytfilt(indici_dirt) = 0;
        output.yant = ifft(output.ytfilt);
        end
end

end
