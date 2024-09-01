function result = SC(s)

% SC          Calculates the resulting value
%             for each beta and starting index
% 
% -- Input --
% * data      Input signal (raw array and not processed)
%
% * N         Number of subcarriers
%
% * Ng        Cyclic prefix length
%
% * Fs        Sampling frequency in the receiver
%
% * bs        Bits per symbol
%
% * FcEst     A priori estimate of Fc and
%             the exact center of the band captured to produce the input
%
%   -- Output -- 
% * result    The value for maximum comparison
%


% The following data should be provided as input:
input = s.data;

N = s.N;
Ng = s.Ng;
Fs = s.Fs;

% 2 for 4QAM & 4 for 16QAM
bs = s.bs;
FcEst = s.FcEst;

% The variable parameters
startIndex = s.startIndex;


beta = s.beta;



ty = (0:(N + Ng - 1))/Fs;

y = input((startIndex):(startIndex + Ng + N - 1));
[yResampled , tyResampled] = resample(y , ty , (1 - beta)*Fs);
expo = 2i*pi*beta*FcEst;


if(length(yResampled) < (N + Ng))
    yResampled = [yResampled ; mean(yResampled)*ones(N + Ng - length(yResampled) , 1)];
    timetemp = tyResampled * (1 - beta)*Fs;
    timetemp = [timetemp , length(timetemp):(N + Ng)];
    tyResampled = timetemp/((1-beta)*Fs);
end

% yfixed = zeros(N + Ng , 1);
% for i = 1:(N + Ng)
%     yfixed(i) = yResampled(i)*exp(expo*tyResampled(i));
% end

yfixed = zeros(N + Ng , 1);
for i = 1:(N + Ng)
    yfixed(i) = yResampled(i)*exp(expo*tyResampled(i));
end


data = y((Ng):(Ng + N -1));
Y = fft(data);


seperatedMat = [real(Y) , imag(Y)];
[idx , C , sumd , D] = kmeans(seperatedMat , 2^bs);
twice = D.^2;
len = size(twice , 1);
variance = sum(twice , 1) / (len - 1);
centeroid = C(: , 1) + 1i*C(: , 2);



temp = zeros(2^bs , 1);
for i = 1:2^bs
    temp(i) = (abs(centeroid(i))^2)/(2*variance(i));
end
result = mean(temp);


end