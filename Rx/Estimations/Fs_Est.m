function Fs = Fs_Est(s)

% N_Est          Estimates the # of subcarriers using a received signal
% 
% -- Input --
% * data      Input signal (raw array and not processed)
%
% * N         The Estimated number of subcarriers
%
% * Fsr       Samplig frequency in the receiver
%
%   -- Output -- 
% * Fs         Estimated Fs parameter
%

% The following data should be provided as input
N = s.N;
input = s.data;

% Turn on if you need to see the data autocorrelateion plot
if(~isfield(s , 'info'))
    s.info = true;
end
info = s.info;

% The guessed Fs based on the power spectrum
if(~isfield(s , 'FsGuess'))
    s.FsGuess = 2.2e8;
end
FsGuess = s.FsGuess;

% The possible error between guessed Fs and actual Fs
% (Shown as p in the paper)
if(~isfield(s , 'err'))
    s.err = 1/10;
end
err = s.err;

% sureValue is the number of elements in each correlation try
% This number should be big enough to contain at least (N + Ng)
% 5316 is the number based on the estimations using Cramer-Rao bound
if(~isfield(s , 'sureValue'))
    s.sureValue = 5316;
end
sureValue = s.sureValue;

eta = s.Fsr/FsGuess;
Srb = ceil(N*eta*(1-err)):floor(N*eta*(1+err));
% In order for reshape to work the data needs to be truncated
inputLength = length(input);
signal = input(1:floor(inputLength/sureValue)*sureValue);

% Reshaping the data to be able to process any amount of data
vecmat = reshape(signal , sureValue , []);
[p , q] = size(vecmat);

vec = zeros(2*p-1 , 1);
for i = 1:q
    this = vecmat(: , i);
    temp = xcorr(this,this);
    vec = vec + temp;
end

% The peak is for zero delay in correlation
[val1 , idx1] = max(abs(vec));

results = abs(vec(idx1 + Srb));
[valNr , idxNr] = max(results);
Nr = Srb(1) + idxNr - 1;
Fs_Estimate = round((N*s.Fsr)/(Nr*1e6))*1e6;

if(info)
    Nr
    Fs_Estimate
end

Fs = Fs_Estimate;
end