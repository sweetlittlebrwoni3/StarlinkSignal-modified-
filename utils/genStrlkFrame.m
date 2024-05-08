function frame = genStrlkFrame(s)
% genStrlkFrame     generates a Starlink frame, with data packaged into frames 
%                   transmitted at 750 Hz consisting of 302 symbols, two
%                   of which are known
% 
% -- Input --
% A struct with the following fields:
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
% * type      (optional) To be provided if *data is also provided. 
%             Modulation type ('PSK' and 'QAM' only, for PSK only M=4 is allowed)
%             So far, QAM and 4QAM with pi/4 offset have been observed 
%             (in Matlab PSK with M = 4). 
%
% * doppVec   (optional) phase history for doppler applied to stream.
%             Mutually exclusive with beta input. Provide the phase history
%             sampled at 240 MHz
%
% * sigma2_w  (optional) Noise variance. Taken with the required SNR input,
%             this determins the noise amplitude and signal amplitude.
%
% * channel   (optional) A comms toolbox channel object.
%
%   -- Output -- 
%
%   y       Starlink OFDM frame in time, sampled at Fsr, centered at Fsr,
%           expressed at baseband.

SNR = 10^(s.SNRdB/10);
sigma2_w = 1/SNR;
A = 1;
if (isfield(s,'sigma2_w'))
    sigma2_w = s.sigma2_w;
    A = sqrt(SNR.*sigma2_w);
end


Tfg = (68/15)*1e-6; % Frame guard interval
s.Fs = 240e6; % Starlink signal BW
Nfg = round(Tfg*s.Fs);
s.N = 1024; % Starlink # of subcarriers
Fdelta = 250e6; % Starlink channel spacing
% Get closest starlink channel center to receiver center
Fcii = getClosestFch(s.Fcr);

SNRdB = s.SNRdB;
beta = s.beta;
Fsr = s.Fsr;
Fcr = s.Fcr;

s.SNRdB = nan();
s.beta = 0;
s.Fcr = Fcii;
s.Fsr = s.Fs;

%-------------------------------------------%
%--------Signal at TX without noise---------%
%-------------------------------------------%
% The starlink frame is made up of 303 subframes.
% The first subframe is the PSS (Primary Sync Sequence), index 0
PSS = genPss();
% The second subframe is the SSS (Secondary Sync Seq), index 1
SSS = genSss();

if (~isfield(s,'data'))
    % Frames with index around 2-5 are 4PSK. Choosing 4.
    s.Nsym = 4;
    s.type = 'PSK';
    s.Midx = 4;
    Data = genStrlkOFDM(s);
    % Frames with index around 6-9 are 16QAM. Choosing 4.
    s.Nsym = 4;
    s.type = 'QAM';
    s.Midx = 16;
    Data = [Data; genStrlkOFDM(s)];
    % Frames with index around 10-300 are 4QAM. Choosing remaining 292
    s.Nsym = 292;
    s.type = 'QAM';
    s.Midx = 4;
    Data = [Data; genStrlkOFDM(s)];    
else
    s.Nsym = 300;
    Data = genStrlkOFDM(s);
end
% The last subframe is empty acting as a frame guard, index 302
Fg = zeros(Nfg,1);

 
s.Fcr = Fcr;
s.Fsr = Fsr;
s.beta = beta;
s.SNRdB = SNRdB;

frame = [PSS; SSS; Data; Fg];

%----- Simulate Doppler & Receiver center bias----%
if (isfield(s,'beta'))
        if (isfield(s,'channel'))
            % channel = s.channel();
            % frame = s.chanFilt(frame,channel);
            FD = -s.beta*getClosestFch(s.Fcr); % Doppler
            release(s.channel)
            s.channel.DirectPathDopplerShift = FD;
            Nsym = length(frame)/1056;
            for ii = 1:Nsym
                iidum = ((ii-1)*1056+1):(ii*1056);
                frame(iidum) = s.channel(frame(iidum));
            end
        else
            FD = -s.beta*getClosestFch(s.Fcr); % Doppler
            fShift = FD + getClosestFch(s.Fcr) - s.Fcr; % total frequency shift

            tVec = (0:length(frame)-1)'/s.Fs;
            frame = frame.*exp(1j*2*pi*fShift*tVec);
        end
elseif (isfield(s,'doppVec'))
        if (isfield(s,'channel'))
            fShift = getClosestFch(s.Fcr) - s.Fcr; % total frequency shift
            offsetVec = 2*pi*(s.doppVec-s.doppVec(1))./Fs;
            tVec = (0:length(y)-1)'/Fs;
            
            release(s.channel)
            s.channel.DirectPathDopplerShift = s.doppVec(1);
            Nsym = length(frame)/1056;
            for ii = 1:Nsym
                iidum = ((ii-1)*1056+1):(ii*1056);
                frame(iidum) = s.channel(frame(iidum));
            end
            % channel = s.channel();
            % frame = s.chanFilt(frame,channel);

            Phihist = cumsum(offsetVec);    
            frame = frame.*exp(1j*(Phihist + 2*pi*fShift*tVec));        
                
        else
            fShift = getClosestFch(s.Fcr) - s.Fcr; % total frequency shift
            tVec = (0:length(y)-1)'/s.Fs;
            offsetVec = 2*pi*s.doppVec./s.Fs;
            Phihist = cumsum(offsetVec);
            frame = frame.*exp(1j*(Phihist + 2*pi*fShift*tVec));
        end
end


frame = A.*frame;

%----- Simulate Noise for PSS, SSS and FG  ----%
if (~isnan(s.SNRdB))    
    sigmaIQ = sqrt(sigma2_w/2);
    Nsamps = length(frame);
    nVec = complex(sigmaIQ*randn(Nsamps,1),sigmaIQ*randn(Nsamps,1));
    tmp = frame;
    frame = frame + nVec;
end

%----- Resample if necessary ----%
if(s.Fsr < s.Fs)
    tVec = (0:length(frame)-1)'/s.Fs;
    frame = resample(frame,tVec,s.Fsr,"spline");
end

end