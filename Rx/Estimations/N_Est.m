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



% The guessed Fs based on the power spectrum
FsGuess = 2.2e8;
% The possible error between guessed Fs and actual Fs
% (Shown as p in the paper)
err = 1/10;
% sureValue is the number of elements in each correlation try
% This number should be big enough to contain at least (N + Ng)
% 5316 is the number based on the estimations using Cramer-Rao bound
sureValue = 5316;

S = 2.^(9:12);
% In order for reshape to work the data needs to be truncated
inputLength = length(s.data);
signal = s.data(1:floor(inputLength/sureValue)*sureValue);

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
N_Estimate = abs(idx2 - idx1)

t = 1:length(vec);
plot(t,mag2db(abs(vec)));

ests = vec(S + idx1);
[tempval , tempidx] = max(ests);
N_newEstimate = S(tempidx)


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
N_bestEstimate = 2^(8 + maxIdx)

N = N_bestEstimate;

end






