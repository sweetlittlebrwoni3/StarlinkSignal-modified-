function [fig] = plotSpec(data,tstart,tdur,Fcr,Fsr,Fc,Fs,NFFT,Stitle)
% plotSpec plots a spectrogram of the provided data. 
%
%   data     Either of vector of complex samples, or a .bin containing
%            I and Q samples, considered the provided signal
%   
%   tstart   Start of plot in seconds relative to beginning of input data
%   
%   tdur     Duration of plot in seconds. If 0 the whole section is
%            plotted only if the input data is a vector of samples.
%   
%   Fcr      Center frequency of Passband equivalent signal (used only
%            if desired plot center is not the same as signal center frq)
%   
%   Fsr      Sampling rate of provided signal
%
%   Fc       Center frequency desired to be plotted (used only
%            if desired plot center is not the same as signal center frq)   
%   
%   Fs       Desired plot sampling rate. Will resample provided signal to
%            Fs if Fs != Fsr
%    
%   NFFT     # of FFT points for spectrogram
%
%   Stitle   Optional title for plot.

    if ~exist('Stitle','var')
     % third parameter does not exist, so default it to something
      Stitle = sprintf("Spectrogram centered at %d",Fc);
    end
   
    if (tdur == 0)
        tdur = length(data)/Fsr;
    end
    seekOffset = floor(tstart*Fsr);
    Nblock = floor(tdur*Fsr);
    data = data(seekOffset+1:seekOffset+Nblock);
    
    tyVec = (0:length(data)-1)'/Fsr;
    % resample to Fs
    if (Fs ~= Fsr)
        [data, tyVec] = resample(data,tyVec,Fs);
    end
    % Shift by RX bias
    if (Fc ~= Fcr)
        Fshift = Fcr - Fc;
        data = data.*exp(1i*2*pi*Fshift*tyVec);
    end
    
    spectrogram(data/sqrt(mean(data.*conj(data))),kaiser(NFFT,0.5),...
              floor(NFFT/2),NFFT,Fs,'psd','centered', ...
              'yaxis', 'minthreshold', -90);
    title(Stitle);
end
