function [beta0 , nm0] = Sym_CF_Est(s)
% Sym_CF_Est                Estimates the first element of a symbol and
%                           carrier frequency offset of that symbol
% 
% -- Input --
% * data                    Input signal (raw array and not processed)
%
% * Fsr                     Sampling frequency in the receiver
%
% * N                       Number of subcarriers
%
% * Ng                      Cyclic prefix length
%
% * Tf                      Frame length
%
% * Fcr                     A priori estimate of Fc and
%                           the exact center of the band captured to produce the input
%
% * bs                      Number of bits per symbol
%
%   -- Output -- 
% * nm0                     Estimated first symbol element
%
% * beta0                   Estimated first symbol doppler factor
%




% The switch to see more details
if(~isfield(s , 'info'))
    s.info = false;
end
info = s.info;

% The following data should be provided as input
N = s.N;
Ng = s.Ng;
Fs = s.Fs;
Tf = s.Tf;
% Nsf = s.Nsf;
Fcr = s.Fcr;
input = s.data;


% Parameters to be determined
if(~isfield(s , 'bs'))
    s.bs = 2;
end
bs = s.bs;

if(~isfield(s , 'symbolIndex'))
    s.symbolIndex = 1;
end
symbolIndex = s.symbolIndex;

if(~isfield(s , 'betamEst'))
    s.betamEst = 2.5e-5;
end
betamEst = s.betamEst;

if(~isfield(s , 'betam'))
    s.betam = 2.5e-5;
end
betam = s.betam;

if(~isfield(s , 'FcEst'))
    s.FcEst = Fcr;
end
FcEst = s.FcEst;


% Expected error in nmi0 approximation
d = 10;

% Approximation using signal energy increase observation
% in the begining of the frame ( y(n)^2 )
data = input(1:2*(Fs*Tf));
% The minus one is only based on the impirical results
nm00Approx = findchangepts(abs(data)) - 1;

if(info)
    findchangepts(abs(data));
end

% Commented till further inspections
% % Set of all the first samples of OFDM symbols in a frame
% nm0Approx = (N + Ng)*(0:Nsf) + nm00Approx;
% 
% Sm = zeros(length(nm0Approx) , 2*d + 1);
% for i = 1:length(nm0Approx)
%     Sm(i , :) = (-d:d) + nm0Approx(i);
% end

% Test scenario for one symbol
nm01Approx = nm00Approx + (N + Ng)*symbolIndex + (-d:d);


s.FcEst = FcEst;
s.data = data;
s.bs = bs;

% Epsilon should be limited to a few percents
epsilon = 0.01;
deltaBeta = (epsilon * Fs)/(N * FcEst);
q = ceil((betamEst - betam)/deltaBeta):floor((betamEst + betam)/deltaBeta);
Bm = q * deltaBeta;


resmat = zeros(length(nm01Approx) , length(Bm));
for i = 1:length(nm01Approx)
    s.startIndex = nm01Approx(i);
    for j = 1:length(Bm)
        s.beta = Bm(j);
        resmat(i , j) = SC(s);
    end
end
this = reshape(resmat' , 1 , []);
[maxVal , maxIdx] = max(this);
nm00 = nm01Approx(floor(maxIdx/length(Bm))) - (N + Ng + 1);
betaidx = rem(maxIdx - 1 , length(Bm)) + 1;
betam = Bm(betaidx);

nm0 = nm00;
beta0 = betam;


end