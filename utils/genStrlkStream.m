function [y,present] = genStrlkStream(s)
% genStrlkStream    generates a stream of Starlink frames, with data packaged 
%                   into frames transmitted at 750 Hz consisting of 302 symbols, 
%                   two of which are known
% -- Input --
% A struct with the following fields:
%
% * Midx      Modulation index (2 = BPSK, 4 = 4QAM/4PSK, 16 = 16QAM/16PSK, etc.)
%             So far, only M = 4 and M = 16 have been observed.
%  
% * type      Modulation type ('PSK' and 'QAM' only, for PSK only M=4 is allowed)
%             So far, QAM and 4QAM with pi/4 offset have been observed 
%             (in Matlab PSK with M = 4)
%  
% * SNRdB     Simulated Signal to noise ratio, in dB. For no noise, use nan()
%  
% * Fsr       Receiver sample rate.  If it is less than Fs, the signal is filtered
%             before resampling to prevent aliasing.  Set Fsr to Fs to skip 
%             resampling and get the full Fs signal.
%   
% * Fcr       Receiver center frequency, in Hz
%  
% * beta      Doppler factor. Doppler shift is FD = -beta*Fc where 
%             beta = vlos/c. Simulated by resampling at (1+beta)*Fs and
%             shifting by FD. Set beta to zero for no Doppler shift.
%
% * Tdur      Duration of stream
%  
% * data      1024 x K vector (optional). Each column corresponds to the
%             serial data transmitted on a symbol's subcarriers. Each column
%             should have elements in the range [0 , M-1], where M determines 
%             the constellation size for the symbol. The data per symbol should be provided 
%             such that the first index corresponds to the ceil(N/2) subcarrier, 
%             the element at index ceil(N/2) corresponds to the 1st subcarrier,
%             and the final index corresponds to the floor(N/2) subcarrier value.
%             The number of symbols whose data can be set K should be in the
%             range [1 , Nsym]. If K < Nsym, the
%             symbols will take on the columns 1024xK and wrap for the
%             remaining 1024x(1:(Nsym-K)) to be the first Nsym-K columns of
%             data.
%
% * prob      (optional) probability of frame being present in a frame
%             slot. If not passed through, the a frame will always be present 
%             (aka consecutive frames)
%
% * present   (optional) Kx1 vector of 1s or 0s indicating a frame present
%             at cadences of 1/750s. Mutually exlusive with prob option. If present
%             option is chosen, Tdur is ignored. 
%
% * tau       (optional) time offser from first frame slot. No frame is
%             placed in partial first frame slot.
%
% * doppVec    (optional) phase history for doppler applied to stream.
%             Mutually exclusive with beta input. Provide the phase history
%             sampled at 240 MHz
%
% * sigma2_w  (optional) Noise variance. Taken with the required SNR input,
%             this determins the noise amplitude and signal amplitude.
%
%   -- Output -- 
%
%   y       Starlink OFDM symbol in time, sampled at Fsr, centered at Fsr,
%           expressed at baseband.

SNR = 10^(s.SNRdB/10);
sigma2_w = 1/SNR;
A = 1;
if (isfield(s,'sigma2_w'))
    sigma2_w = s.sigma2_w;
    A = sqrt(SNR.*sigma2_w);
end

if (isfield(s,'present'))
    if (isfield(s,'prob'))
        error("Can specify either prob or present, not both.");
    end
end
if (~isfield(s,'prob'))
    prob = 1;
else 
    prob = s.prob;
end
if (~isfield(s,'tau'))
    tau = 0;
else 
    tau = s.tau;
end

Fs = 240e6; % Starlink signal BW
Tframe = 1/750; % Starlink Frame duration
Nfr = ceil(s.Tdur/Tframe); % # of whole Starlink Frames in stream
if (isfield(s,'present'))
    present = s.present;
    Nfr_p = length(s.present);
    if (Nfr>Nfr_p)
        present(end+1:end+Nfr-Nfr_p) = 0;
    else
        present = present(1:Nfr);
    end
else
    % Indication whether starlink frame is present in slot
    present = binornd(1,prob,1,Nfr);
end


Ns = floor(s.Tdur*Fs); % # of samples in stream
Nframe = floor(Tframe*Fs); % # samples per full BW Starlink frame

if (isfield(s,'beta'))
    if (isfield(s,"phHist"))
        error("You can only specify a CFO (beta) OR a phase time-history, not both.")
    end
end
if (isfield(s,"phHist"))
    if (isfield(s,'beta'))
        error("You can only specify a CFO (beta) OR a phase time-history, not both.")
    end
    if (length(s.phHist) ~= Ns)
        error("Your phase time-history should be %d long according to your inputs.",Ns)
    end
end
    
y = zeros(Ns,1);



for ii = 1:Nfr
    if (present(ii)) % Place frame
        fr.SNRdB = nan(); % Signal to noise ratio, in dB
        fr.Fsr = Fs; % Receiver sample rate, in Hz. 
        fr.Fcr = getClosestFch(s.Fsr); % Receiver center frequency, in Hz
        fr.beta = 0; % Doppler parameter
        if (isfield(s,'type'))
             fr.type = s.type; % subcarrier constellation type
        end
        if (isfield(s,'data'))
             fr.data = s.data; % subcarrier constellation type
        end
        frame = genStrlkFrame(fr);

        y((ii-1)*Nframe+1:((ii-1)*Nframe+length(frame))) = frame;
    end
end

Ntau = round(s.tau.*Fs);
y = [zeros(Ntau,1);y];

y = A.*y;

%----- Simulate Doppler & Receiver center bias----%
if (isfield(s,'beta'))
    if(~(s.beta == 0 && getClosestFch(s.Fcr) == s.Fcr))
        FD = -s.beta*getClosestFch(s.Fcr); % Doppler
        fShift = FD + getClosestFch(s.Fcr) - s.Fcr; % total frequency shift
    
        tVec = (0:length(y)-1)'/Fs;
        y = y.*exp(1j*2*pi*fShift*tVec);
    end
elseif (isfield(s,'doppVec'))
            fShift = getClosestFch(s.Fcr) - s.Fcr; % total frequency shift
            tVec = (0:length(y)-1)'/Fs;
            offsetVec = 2*pi*s.doppVec./Fs;
            Phihist = cumsum(offsetVec);
            y = y.*exp(1j*(Phihist + 2*pi*fShift*tVec));
end


%----- Simulate Noise for PSS, SSS and FG  ----%
if (~isnan(s.SNRdB))    
    sigmaIQ = sqrt(sigma2_w/2);

    Nsamps = length(y);
    nVec = complex(sigmaIQ*randn(Nsamps,1),sigmaIQ*randn(Nsamps,1));
    y = y + nVec;
end

%----- Resample if necessary ----%
if(s.Fsr < Fs)
    tVec = (0:length(y)-1)'/Fs;
    y = resample(y,tVec,s.Fsr,"spline");
end


end