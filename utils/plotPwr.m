function [fig] = plotPwr(data,Nblocks,tstart,tdur,Fs,NFFT)
%   plotPwr plots the PSD and Power vs Time of the provided signal
%           starting at tstart for a duration of tdur.
%   -- Input --
%
%   data     vector of complex samples
%
%   Nblocks  # of blocks to sweep for power vs time plot
%
%   tstart   Start of data to process relative to the first sample
%
%   tdur     Duration of data to process starting from to tstart
%   
%   Fs       Sampling rate of provided signal
%
%   NFFT     # of FFT points for PSD plot
%

    if (Nblocks>10*tdur*Fs)
        error("Provided number of blocks is too high for duration specified.")
    end
    seekOffset = floor(tstart*Fs);
    Nblock = floor(tdur*Fs);
    data = data(seekOffset+1:seekOffset+Nblock);
    

    tVec = (0:(length(data)-1))'/Fs+tstart;
    [Syy,fVec] = pwelch(data,kaiser(NFFT,3),NFFT/2,NFFT,Fs,'psd','centered');
    subplot(2,1,1)
    plot(fVec/1e6,10*log10(Syy));
    grid on;
    xlabel('Frequency (MHz)');
    ylabel('Power density (dB/Hz)');
    title('Power spectral density estimate');
    L_seg = floor(length(data)/Nblocks);
    Pwrdb = zeros(Nblocks,1);
    for i=1:Nblocks
        [Syy,~] = pwelch(data(L_seg*(i-1)+1:L_seg*(i)),kaiser(round(L_seg/2),3),[],NFFT,Fs,'psd','centered');
        Pwrdb(i) = max(10*log10(sum(Syy.*Fs)),-180); % dB/Hz * Hz = dB
    end
    subplot(2,1,2)
    tVec = linspace(0,tdur,Nblocks)+tstart;
    plot(tVec,Pwrdb);
    grid on;
    xlabel('Time (s)');
    ylabel('dB');
    title('Power vs Time');
end
