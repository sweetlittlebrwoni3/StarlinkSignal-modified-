function y = genStrlkOFDM(s)
% genStrlkOFDM is an OFDM symbol generator with specific Starlink parameters. 
% 
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
%  * Nsym      (optional) Number of consequtive symbols to generate.
%              Default is 1
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
%   -- Output -- 
%
%   y       Starlink OFDM symbol in time, sampled at Fsr, centered at Fsr,
%           expressed at baseband.

% Channel bandwidth (not including channel guard intervals), symbol rate of
% the serial data sequence, in Hz.
s.Fs = 240e6;
% Number of subcarriers in bandwidth Fs, size of FFT/IFFT, number of
% samples at rate Fs in the useful symbol interval. For OFDM, must be a
% power of 2.
s.N = 1024;
% Duration of OFDM symbol guard interval (cyclic prefix), expressed in
% intervals of 1/Fs.  Typically a power of 2.  Never more than half of N.
s.Ng = 32;
% Enable 4 subcarrier gutter at center
s.gutter = 1;

if (~isfield(s,'Nsym'))
    s.Nsym = 1;
end
if (isfield(s,'data') && ~isfield(s,'type'))
    error(" A constellation type should also be provided as s.type (i.e. 'PSK' or 'QAM')")
end
if (isfield(s,'data'))
[l,w] = size(s.data);
if(l ~= s.N || w>300)
    error('s.data must be an N x 300 vector at most');
end
end



if((~isfield(s,'data')) && (mod(log2(s.Midx),2) ~= 0) && (s.Midx ~= 2))
    error('Midx must 2 or an even power of 2');
end

% Get closest starlink channel center to receiver center
F = s.Fs/s.N;
chIdx = round((s.Fcr/1e9 - 10.7 - F/2/1e9)/0.25 + 0.5);
Fcii = (10.7e9 + F/2 + 250e6*(chIdx - 0.5));
s.Fc = Fcii;

y = genOFDM(s);

end
