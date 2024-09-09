function Tf = Tf_Est(s)

% N_Est          Estimates the # of subcarriers using a received signal
% 
% -- Input --
% * data      Input signal (raw array and not processed)
%
% * Fs        Signal bandwidth
%
% * N         Number of subcarriers
%
% * Ng        Cyclic prefix length
%
%   -- Output -- 
% * Tf         Estimated Tf parameter
%

% ** There needs to be at least two frames together for this estimation to
% work

% The following data should be provided as input
Fs = s.Fs;
N = s.N;
Ng = s.Ng;
input = s.data;


% The switch to see more details
if(~isfield(s , 'info'))
    s.info = false;
end
info = s.info;


% The obtained upper bound on the
% smallest active singal interval
if(~isfield(s , 'Tm'))
    s.Tm = length(input)/Fs;
end
Tm = s.Tm;


% This field decides between the original
% and modified versions of this Estimation
if(~isfield(s , 'modified'))
    s.modified = true;
end
modified = s.modified;


temp = xcorr(input , input);


% This modification is made to make sure
% the first maxima is chosen
if(modified)
    for i = 1:length(temp)
        temp(i) = temp(i)/i;
    end
end


[maxVal , maxIdx] = max(temp);
Sf = (N + Ng + 1):floor(Fs*Tm - 1);
vec = temp(Sf + maxIdx);
[maxVal2 , maxIdx2] = max(abs(vec));
Nf = Sf(maxIdx2);
result = 1/round(Fs/Nf);

if(info)
    Nf
    result
end

Tf = result;

end