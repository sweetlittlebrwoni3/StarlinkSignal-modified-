function N = N_Est(s)

% N_Est          Estimates the # of subcarriers using a received signal
% 
% -- Input --
% * data      Input signal (raw array and not processed)
%
% * Fsr       Sampling frequency in the receiver
%
%   -- Output -- 
% * N         Estimated N parameter
%

% The following data should be provided as input
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

% The possible range for N
S = 2.^(9:12);

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
vec(idx1) = 0;
[val2 , idx2] = max(abs(vec));
N_Estimate = abs(idx2 - idx1);
if(info)
    N_Estimate
end

% Plotting the autocorrelation in order to see the peaks
if(info)
    t = 1:length(vec);
    plot(t , mag2db(abs(vec)));
end

ests = vec(S + idx1);
[tempval , tempidx] = max(ests);
N_newEstimate = S(tempidx);
if(info)
    N_newEstimate
end


% Creating the Sr dataset
eta = s.Fsr/FsGuess;
maxDim = ceil(max(S)*eta*2*err);
Sr = zeros(length(S) , maxDim);
for i = 1:length(S)
    Srb = ceil(S(i)*eta*(1-err)):...
        floor(S(i)*eta*(1+err));
    for j = 1:length(Srb)
        Sr(i,j) = Srb(j);
    end
end

% The Validation Test
results = abs(vec(Sr + idx1));
Mx = max(results , [] , 2);
valid = false;
while(~valid)
    [maxVal , maxIdx] = max(Mx);
    finishLine = round(S(maxIdx)*eta*2*err);
    [minVal , minIdx] = min(results(maxIdx , 1:finishLine));
    if(mag2db(maxVal) - mag2db(minVal) > 10)
        valid = true;
    elseif(~isempty(Mx))
        Mx(MaxIdx) = [];
    else
        error('Could not estimate N. Please try again.');
    end
end

N_bestEstimate = 2^(8 + maxIdx);
if(info)
    N_bestEstimate
end

N = N_bestEstimate;

end






