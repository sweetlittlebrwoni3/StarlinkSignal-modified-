function N = N_Est(signal)

% N_Est          Estimates the # of subcarriers using a received signal
% 
% -- Input --
% * data      Input signal (raw array and not processed)
%
%   -- Output -- 
% * N         Estimated N parameter
%

dOFDM = (bs*Fs*N)/(N+Ng);

Tg = Ng/Fs;

F = Fs/N;

q = 9:12;
S = 2.^q;


end






