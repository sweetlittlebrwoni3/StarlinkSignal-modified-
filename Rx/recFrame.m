function data = recFrame(s)

% recFrame          recovers the data from a Starlink frame
% 
% -- Input --
% A struct with the following fields:
%  
% * Fsr       Receiver sample rate.  If it is less than Fs, the signal is filtered
%             before resampling to prevent aliasing.  Set Fsr to Fs to skip 
%             resampling and get the full Fs signal.
%   
% * Fcr       Receiver center frequency, in Hz
%
% * y         Starlink OFDM frame in time, sampled at Fsr, centered at Fsr,
%             expressed at baseband.
%
% * type      (optional) To be provided if *data is also provided. 
%             Modulation type ('PSK' and 'QAM' only, for PSK only M=4 is allowed)
%             So far, QAM and 4QAM with pi/4 offset have been observed 
%             (in Matlab PSK with M = 4).
%
% * channel   (optional) A comms toolbox channel object.
%
%   -- Output -- 
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



Tfg = (68/15)*1e-6; % Frame guard interval
s.Fs = 240e6; % Starlink signal BW
Nfg = round(Tfg*s.Fs);
pss = genPss();
sss = genSss();
input = s.y(length(pss)+length(sss)+1:length(s.y)-Nfg);


ss.type = s.type;
ss.Midx = s.Midx;
ss.Fsr = s.Fsr;
ss.Fs = s.Fs;
ss.N = s.N;
ss.Ng = s.Ng;
ss.Fc = s.Fc;
ss.Fcr = s.Fcr;
ss.gutter = s.gutter;
ss.Nsym = s.Nsym;
ss.beta = s.beta;

input1 = reshape(input,s.N+s.Ng,[]);

data = zeros(s.N - 4,s.Nsfd);
for ii = 1:s.Nsfd
    if(~isfield(s,"type"))
        if(ii < 5)
            ss.type = "PSK";
            ss.Midx = 4;
        else
            ss.type = "QAM";
            ss.Midx = 16;
        end
    end
    ss.yVec = input1(:,ii);
    sym = recOFDM(ss);
    data(:,ii) = sym;
end






