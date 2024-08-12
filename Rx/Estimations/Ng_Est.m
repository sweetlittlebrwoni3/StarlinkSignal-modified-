function Ng = Ng_Est(s)

% N_Est          Estimates the # of subcarriers using a received signal
% 
% -- Input --
% * data      Input signal (raw array and not processed)
%
% * N         The Estimated number of subcarriers
%
% * Fs        The Estimated Bandwidth
%
%   -- Output -- 
% * Ng         Estimated Ng parameter
%


% The following data should be provided as input
input = s.data;
N = s.N;
Fs = s.Fs;


% The option to choose amongst the methods
if(~isfield(s , 'Ng_method'))
    s.Ng_method = "optimized";
end
method = s.Ng_method;


% Worst-case 95% root-mean-square delay spread
% for the Ku-band
if(~isfield(s , 'Td'))
    s.Td = 108e-9;
end
Td = s.Td;


% No improvement attains to values of Np above N/Ng
% Np = ceil(N/min(2*q)); can be used to be sure
% But if the number of samples is bigger than the 
% number of samples in a frame, 1 would be fine for Np
if(~isfield(s , 'Np'))
    s.Np = 1;
end
Np = s.Np;

% This works best if M > Tf*Fs and SNR > 3.5 dB
if(~isfield(s , 'processLength'))
    s.processLength = length(input);
end
M = s.processLength;


% Constructing the Sg set to select (Ng + N) from
b = 2.^(9:12);
% Constructing the Sg dataset
q = ceil(Td*Fs/4):floor(Td*Fs);
Sg = zeros(length(q) , length(b));
for i = 1:length(q)
    for j = 1:length(b)
        Sg(i , j) = 2*q(i) + b(j);
    end
end
zeta = reshape(Sg,1,[]);
% Sg dataset is constructed
% but since the only possible value for b is 1024
% the other values are unused.



shiftedInput = circshift(input , N);


% This method is faster but has less accuracy
if(method == "fast")
    zeta = 1024 + 2*q;
    vec = zeros(length(input) , 1);
    for i = 1:length(input)
        vec(i) = shiftedInput(i) * conj(input(i));
    end
    test = xcorr(vec , vec);
    [val , idx] = max(test);
    here = zeros(length(zeta) , 1);
    for i = 1:length(zeta)
        here(i) = test(idx + zeta(i));
    end
    [maxVal , maxIdx] = max(abs(here));
    Ng_fast = zeta(maxIdx) - N;

    Ng = Ng_fast;
end


% This is the method introduced in the paper
if(method == "basic")
    % Since N is estimated, b can be considered to be 1024
    % b = 1024;
    % zeta = 2*q + b;
    alpha = (0:Np)./zeta';

    temp = 0;
    for n = 1:M
        temp = temp +...
            2*shiftedInput(n)*conj(input(n))...
            *cos(-2*pi*(n-1)*alpha);
    end

    sigma2 = sum((1/M)*abs(temp) , 2);

    [MaxValOpt , MaxIdxOpt] = max(sigma2);
    Ng_basic = zeta(MaxIdxOpt) - N;

    Ng = Ng_basic;
end


% This method is optimized to run faster in matlab
if(method == "optimized")
    % Since N is estimated, b can be considered to be 1024
    % b = 1024;
    % zeta = 2*q + b;
    alpha = (0:Np)./zeta';

    n = 0:(M-1);
    n_extended = reshape(n, [1, 1, numel(n)]);
    this = alpha .* n_extended;
    freq = cos(-2*pi*this);
    hit = 2*shiftedInput.*conj(input);
    temp = reshape(hit, [1 , 1 ,numel(input)]);
    result = freq .* temp;
    sigma = sum((1/M)*abs(result) , 3);
    sigma2 = sum(sigma , 2);

    [MaxVal2 , MaxIdx2] = max(sigma2);
    Ng_optimized = zeta(MaxIdx2) - N;

    Ng = Ng_optimized;
end


end