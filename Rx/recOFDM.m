function data = recOFDM(s)
% recOFDM recovers the data from an OFDM symbol
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
%  * yVec      Input OFDM symbol
%              
%   -- Output -- 
%
%   data       Demodulated OFDM symbol back to it's normal data type
%


% Checking to see if there's gutter
if(~isfield(s,'gutter'))
    Ngut = 0;
else
    if(s.gutter == 0)
        Ngut = 0;
    else
    Ngut = 4;
    end
end
Nd = s.N - Ngut;

% Removing the cyclic prefix
input = s.yVec(s.Ng+1:length(s.yVec));


%input = fftshift(input);

% Demodulation using ofdmdemod:
x1 = ofdmdemod(input,1024,0);
x1 = x1(Ngut/2+1:s.N-Ngut/2);
x1 = fftshift(x1);

if ( upper(string(s.type)) == "PSK")
    x1Vec = pskdemod(x1,s.Midx);
elseif( upper(string(s.type)) == "QAM")
    x1Vec = qamdemod(x1,s.Midx,'UnitAveragePower',true);
else
    error('Type must be "QAM" or "PSK".');
end


% Demodulation using fft:
x2 = (1/sqrt(s.N))*fft(input);
x2 = fftshift(x2);
x2 = x2(Ngut/2+1:s.N-Ngut/2);
x2 = fftshift(x2);

if ( upper(string(s.type)) == "PSK")
    x2Vec = pskdemod(x2,s.Midx);
elseif( upper(string(s.type)) == "QAM")
    x2Vec = qamdemod(x2,s.Midx,'UnitAveragePower',true);
else
    error('Type must be "QAM" or "PSK".');
end


% Choosing between the two methods
xVec = x1Vec;
if(isfield(s,'method'))
    if(s.method == "fft")
        xVec = x2Vec;
    end
end


% Data recovery
output = zeros(1,Nd);
for ii=1:Nd
    if(real(xVec(ii)) >= 0)
        output(ii) = output(ii) + 1;
    else
        output(ii) = output(ii) - 1;
    end
    if(imag(xVec(ii)) >= 0)
        output(ii) = output(ii) + 1i;
    else
        output(ii) = output(ii) - 1i;
    end
end

data = ifftshift(xVec);