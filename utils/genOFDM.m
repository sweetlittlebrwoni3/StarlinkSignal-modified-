function yVec = genOFDM(s)
% genOFDM is a OFDM signal generator. 
%
% -- Input --
%
% A struct with the following fields:
%
%  * Fs        Channel bandwidth (not including channel guard intervals) in Hz
%
%  * N         Number of subcarriers, size of FFT/IFFT w/out CP (cyclic prefix)
%
%  * Ng        Number of guard subcarriers CP (cyclic prefix)
% 
%  * Midx      Subcarrier constellation size (2 = BPSK, 4 = 4QAM, 16 =
%              16QAM, etc.). OVERWRITTEN IF DATA IS PROVIDED see *data below
%
%  * type      Modulation type ('PSK' or 'QAM' only)
%
%  * SNRdB     Simulated Signal to noise ratio, in dB. For no noise, pass
%              in nan()
%
%  * Fsr       Receiver sample rate.  If it is less than Fs, the signal is filtered
%              before resampling to prevent aliasing.  Set Fsr to Fs to skip 
%              resampling and get the full Fs signal.
%
%  * Fcr       Receiver center frequency.
% 
%  * Fc        OFDM signal center frequency, in Hz
%
%  * beta      Doppler factor. Doppler shift is FD = -beta*Fc where 
%              beta = vlos/c. Simulated by resampling at (1+beta)*Fs and
%              shifting by FD. Set beta to zero for no Doppler shift. For 
%              approaching SVs (beta < 0), the measured inter-frame interval 
%              is shorter (compressed) compared to the ideal inter-frame 
%              interval of 1/750. Assuming a polynomial model of the form
%              tF(t) = p3*t^2 + p2*t + p1, then dtF/dt|_0 = p1 = beta(0).  
%              Doppler frequency shift for a tone at Fc is given by FD = -beta*Fc.
%
%  * gutter    boolean 1 or 0 (optional). Enables a gutter of 4F at center 
%              as observed in Starlink signals
%
%  * Nsym      (optional) Number of consequtive symbols to generate.
%              Default is 1, or if *data is provided, its length.
%
%  * data      1024 x K vector (optional). Each column corresponds to the
%              serial data transmitted on a symbol's subcarriers. Each column
%              should have elements in the range [0 , M-1], where M determines 
%              the constellation size for the symbol. The data per symbol should be provided 
%              such that the first index corresponds to the ceil(N/2) subcarrier, 
%              the element at index ceil(N/2) corresponds to the 1st subcarrier,
%              and the final index corresponds to the floor(N/2) subcarrier value.
%              The number of symbols whose data can be set K should be in the
%              range [1 , Nsym]. If K < Nsym, the
%              symbols will take on the columns 1024xK and wrap for the
%              remaining 1024x(1:(Nsym-K)) to be the first Nsym-K columns of
%              data.
%              
%   -- Output -- 
%
%   yVec       OFDM symbol in time, sampled at Fsr, centered at Fsr,
%              expressed at baseband.
%



%----- Optional parameters & Checks -----------%
if (~isfield(s,'gutter'))
    s.gutter = 0;
end
if (~isfield(s,'Nsym'))
    s.Nsym = 1;
end
% If input value is not even power of 2, then invalid
if(~isfield(s,'data') && (mod(log2(s.Midx),2) ~= 0) && (s.Midx ~= 2))
  error('Midx must 2 or an even power of 2');
end
if (~isfield(s,'data'))
    x = randi([0 s.Midx-1],s.N,s.Nsym);
    [l,w] = size(x);
    s.Midx = 2.^nextpow2(s.Midx).*ones(w,1);
else
    [l,w] = size(s.data);
    if(l ~= s.N || w>s.Nsym)
        error('s.data must be an N x Nsym vector');
    end
    s.Midx = 2.^nextpow2(max(s.data));
    x = s.data;
end
% Generate data symbols each with Midx bits of information


%----- Dependent parameters -----------%
T = s.N/s.Fs; % symbol duration, non cyclic
Tg = s.Ng/s.Fs; % guard duration
Tsym = T + Tg; % ofdm symbol duration
F = s.Fs/s.N; % Subcarrier spacing

%----- Generate simulated serial data symbols -----------%
% Generate those symbols as complex numbers with unit average energy
XVec = zeros(size(x));
if ( upper(string(s.type)) == "PSK")
    for ii = 1:w
        XVec(:,ii) = pskmod(x(:,ii),s.Midx(ii));
    end
elseif( upper(string(s.type)) == "QAM")
    for ii = 1:w
        XVec(:,ii) = qammod(x(:,ii),s.Midx(ii),'UnitAveragePower',true);
    end
else
    error('Type must be "QAM" or "PSK".');
end
if (s.Nsym - w > 0)
    Nreps = floor(s.Nsym/w);
    Nremain = mod(s.Nsym,w);
    XVec = repmat(XVec,[1,Nreps]);
    XVec = [XVec , XVec(:,1:Nremain)];
end

%----- Generate the OFDM symbols from Serial Data ----%
% Organize the complex numbers in batches of N
if (s.gutter)
    Ngut = 4; % Startlink has a 4 subcarrier gutter at center.
    XVec = fftshift(XVec,1);
    XVec((s.N-Ngut)/2+1:(s.N-Ngut)/2+Ngut,:) = 0;
    XVec = fftshift(XVec,1);
end
% Transform to time domain.  Multiply by sqrt(N) to preserve energy
Mx = sqrt(s.N)*ifft(XVec);
% Prepend each symbol with cyclic prefix
MxCP = [Mx(end - s.Ng + 1:end,:); Mx];
% Unroll into serial sample stream
xVec = MxCP(:);

%----- Simulate Doppler & Receiver bias to center----%
if(s.beta ~= 0 || s.Fc ~= s.Fcr)
  tVec = [0:length(xVec)-1]'/s.Fs;
  FD = -s.beta*s.Fc;
  Fshift = FD + s.Fc - s.Fcr;
  xVec = xVec.*exp(j*2*pi*Fshift*tVec);  
end

%----- Pass through AWGN channel
% Since the signal has unit power (unit average energy per symbol), the SNR is
% given by SNR = 1/E[|n(k)|^2] = 1/(2*sigmaIQ^2), where n(k) = nI(k) + j*nQ(k)
% is complex Gaussian noise, with nI and nQ independent zero-mean Gaussian
% noise processes, each having variance sigmaIQ^2.
yVec = xVec;
if (~isnan(s.SNRdB))
    SNR = 10^(s.SNRdB/10);
    sigmaIQ = sqrt(1/(2*SNR));
    Nsamps = length(xVec);
    nVec = complex(sigmaIQ*randn(Nsamps,1),sigmaIQ*randn(Nsamps,1));
    yVec = xVec + nVec;
end

%----- Resample to simulate receiver capture ----%
if(s.Fsr ~= s.Fs && ~isempty(yVec))
  tVec = [0:length(yVec)-1]'/s.Fs;
  yVec = resample(yVec,tVec,s.Fsr);
end
end
