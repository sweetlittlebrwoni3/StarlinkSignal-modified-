clc;
addpath('../Utils')

%--------------------------------------%
%---------OFDM SIGNAL TEST-------------%
%--------------------------------------%
clear all;
% ----------- Test 1 
% No noise, no gutter, Full BW, 1 symbol

N = 1024; % # of s.c
Ng = 32; % c.p length
Midx = 4; % subcarrier constellation size
Nsym = 1; % # of symbols to simulate
Nprovided = 1; % # of provided data sequences to generate


% Test input params
s.Fs = 240e6; % Symbol BW
s.N = N; % # of s.c.
s.Ng = Ng; % c.p length
s.Midx = Midx; % subcarrier constellation size
s.type = 'QAM'; % subcarrier constellation type
s.SNRdB = nan(); % signal to noise ratio, in dB
s.Fsr = s.Fs; % simulated receiver sampling rate
s.Fcr = 11.805e9; % simulated receiver center frequency
s.Fc = 11.805e9; % OFDB Symbol center frequency
s.beta = 0; % CFO
s.gutter = 0; % boolean for gutters
s.data = randi([0 s.Midx-1],s.N,Nprovided); % Input data


% Generate frame, zero pad, and process
y = genOFDM(s);
y = [y; zeros((N+Ng)*ceil(length(y)/(N+Ng)) - length(y),1)]; % zero pad
y = reshape(y,N+Ng,[]); % break up in symbols
y = y(Ng+1:end,:); % remove CP
Y = 1/sqrt(N).*fft(y);

% set gutter values to nan() for easier comparison
if (s.gutter)
    Y(1:2,:) = nan();
    Y(end-1:end,:) = nan();
end

% Expected decoded input data
if ( upper(string(s.type)) == "PSK")
    for ii=1:Nprovided
     seq = pskmod(s.data(:,ii),2^nextpow2(max(s.data(:,ii))));
     if (s.gutter)
     seq(1:2) = nan();
     seq(end-1:end) = nan();
     end
     data_in(:,ii) = getDiffEnc(seq);
    end
elseif( upper(string(s.type)) == "QAM")
    for ii=1:Nprovided
     seq = qammod(s.data(:,ii),2^nextpow2(max(s.data(:,ii))),'UnitAveragePower',true);
     if (s.gutter)
     seq(1:2) = nan();
     seq(end-1:end) = nan();
     end
     data_in(:,ii) = getDiffEnc(seq);
    end
end
if (Nsym - Nprovided > 0)
    Nreps = floor(Nsym/Nprovided);
    Nremain = mod(Nsym,Nprovided);
    data_in = repmat(data_in,[1,Nreps]);
    data_in = [data_in , data_in(:,1:Nremain)];
end

% decoded output data
data_out = getDiffEnc(Y);

% Compare
comp = data_out==data_in;
for ii = 1:Nsym
    match(ii) = sum(isnan(data_in(isnan(data_out(:,ii)),ii)));
end
match = sum(comp) + match;
if (Nsym == sum(match == N))
    str = 'pass';
else
    str = 'Failed';
    warning("Test failed, %d / %d symbols correctly decoded",match,N)
end
fprintf('1. OFDM Test no noise, no gutter, full BW, 1 symbol : %18s \n',str)

clear all;
% ----------- Test 2 
% No noise, no gutter, full BW, Nsym symbols

N = 1024; % # of s.c
Ng = 32; % c.p length
Midx = 4; % subcarrier constellation size
Nsym = 300; % # of symbols to simulate
Nprovided = 57; % # of provided data sequences to generate


% Test input params
s.Fs = 240e6; % Symbol BW
s.N = N; % # of s.c.
s.Ng = Ng; % c.p length
s.Midx = Midx; % subcarrier constellation size
s.Nsym = Nsym; % # of consequtive frames to generate
s.type = 'QAM'; % subcarrier constellation type
s.SNRdB = nan(); % signal to noise ratio, in dB
s.Fsr = s.Fs; % simulated receiver sampling rate
s.Fcr = 11.805e9; % simulated receiver center frequency
s.Fc = 11.805e9; % OFDB Symbol center frequency
s.beta = 0; % CFO
s.gutter = 0; % boolean for gutters
s.data = randi([0 s.Midx-1],s.N,Nprovided); % Input data


% Generate frame, zero pad, and process
y = genOFDM(s);
y = [y; zeros((N+Ng)*ceil(length(y)/(N+Ng)) - length(y),1)]; % zero pad
y = reshape(y,N+Ng,[]); % break up in symbols
y = y(Ng+1:end,:); % remove CP
Y = 1/sqrt(N).*fft(y);

% set gutter values to nan() for easier comparison
if (s.gutter)
    Y(1:2,:) = nan();
    Y(end-1:end,:) = nan();
end

% Expected decoded input data
if ( upper(string(s.type)) == "PSK")
    for ii=1:Nprovided
     seq = pskmod(s.data(:,ii),2^nextpow2(max(s.data(:,ii))));
     if (s.gutter)
     seq(1:2) = nan();
     seq(end-1:end) = nan();
     end
     data_in(:,ii) = getDiffEnc(seq);
    end
elseif( upper(string(s.type)) == "QAM")
    for ii=1:Nprovided
     seq = qammod(s.data(:,ii),2^nextpow2(max(s.data(:,ii))),'UnitAveragePower',true);
     if (s.gutter)
     seq(1:2) = nan();
     seq(end-1:end) = nan();
     end
     data_in(:,ii) = getDiffEnc(seq);
    end
end
if (Nsym - Nprovided > 0)
    Nreps = floor(Nsym/Nprovided);
    Nremain = mod(Nsym,Nprovided);
    data_in = repmat(data_in,[1,Nreps]);
    data_in = [data_in , data_in(:,1:Nremain)];
end

% decoded output data
data_out = getDiffEnc(Y);

% Compare
comp = data_out==data_in;
for ii = 1:Nsym
    match(ii) = sum(isnan(data_in(isnan(data_out(:,ii)),ii)));
end
match = sum(comp) + match;
if (Nsym == sum(match == N))
    str = 'pass';
else
    str = 'Failed';
    warning("Test failed, %d / %d symbols correctly decoded",sum(match == N),Nsym)
end
fprintf('2. OFDM Test no noise, no gutter, full BW, %d symbols : %15s \n',Nsym,str)

clear all;
% ----------- Test 3 
% No noise, with gutter, full BW, Nsym symbols

N = 1024; % # of s.c
Ng = 32; % c.p length
Midx = 4; % subcarrier constellation size
Nsym = 1; % # of symbols to simulate
Nprovided = 1; % # of provided data sequences to generate
epsilon = 1e-6;

% Test input params
s.Fs = 240e6; % Symbol BW
s.N = N; % # of s.c.
s.Ng = Ng; % c.p length
s.Midx = Midx; % subcarrier constellation size
s.Nsym = Nsym; % # of consequtive frames to generate
s.type = 'QAM'; % subcarrier constellation type
s.SNRdB = nan(); % signal to noise ratio, in dB
s.Fsr = s.Fs; % simulated receiver sampling rate
s.Fcr = 11.805e9; % simulated receiver center frequency
s.Fc = 11.805e9; % OFDB Symbol center frequency
s.beta = 0; % CFO
s.gutter = 1; % boolean for gutters
s.data = randi([0 s.Midx-1],s.N,Nprovided); % Input data


% Generate frame, zero pad, and process
y = genOFDM(s);
y = [y; zeros((N+Ng)*ceil(length(y)/(N+Ng)) - length(y),1)]; % zero pad
y = reshape(y,N+Ng,[]); % break up in symbols
y = y(Ng+1:end,:); % remove CP
Y = 1/sqrt(N).*fft(y);

gutter_check = 0;
if (abs(Y(1)) < epsilon && abs(Y(2)) < epsilon && ...
        abs(Y(end)) < epsilon && abs(Y(end-1)) < epsilon)
    gutter_check = 1;
end

% set gutter values to nan() for easier comparison
if (s.gutter)
    Y(1:2,:) = nan();
    Y(end-1:end,:) = nan();
end

% Expected decoded input data
if ( upper(string(s.type)) == "PSK")
    for ii=1:Nprovided
     seq = pskmod(s.data(:,ii),2^nextpow2(max(s.data(:,ii))));
     if (s.gutter)
     seq(1:2) = nan();
     seq(end-1:end) = nan();
     end
     data_in(:,ii) = getDiffEnc(seq);
    end
elseif( upper(string(s.type)) == "QAM")
    for ii=1:Nprovided
     seq = qammod(s.data(:,ii),2^nextpow2(max(s.data(:,ii))),'UnitAveragePower',true);
     if (s.gutter)
     seq(1:2) = nan();
     seq(end-1:end) = nan();
     end
     data_in(:,ii) = getDiffEnc(seq);
    end
end
if (Nsym - Nprovided > 0)
    Nreps = floor(Nsym/Nprovided);
    Nremain = mod(Nsym,Nprovided);
    data_in = repmat(data_in,[1,Nreps]);
    data_in = [data_in , data_in(:,1:Nremain)];
end

% decoded output data
data_out = getDiffEnc(Y);

% Compare
comp = data_out==data_in;
for ii = 1:Nsym
    match(ii) = sum(isnan(data_in(isnan(data_out(:,ii)),ii)));
end
match = sum(comp) + match;
if (Nsym == sum(match == N))
    str = 'pass';
else
    str = 'Failed';
    warning("Test failed, %d / %d symbols correctly decoded",sum(match == N),Nsym)
end
if (~gutter_check)
    str = 'Failed';
    warning("Test failed, gutter was not generated in 4 center subcarriers.")
end
fprintf('3. OFDM Test no noise, with gutter, full BW : %26s \n',str)

clear all;
% ----------- Test 4 
% Noise, with gutter, full BW, Nsym symbols

N = 1024; % # of s.c
Ng = 32; % c.p length
Midx = 4; % subcarrier constellation size
Nsym = 300; % # of symbols to simulate
Nprovided = 1; % # of provided data sequences to generate
epsilon = 1e-6;
SNRdb = 5;
avgSCdecodethresh = 800;
snrthresh_dB = 2;
% Test input params
s.Fs = 240e6; % Symbol BW
s.N = N; % # of s.c.
s.Ng = Ng; % c.p length
s.Midx = Midx; % subcarrier constellation size
s.Nsym = Nsym; % # of consequtive frames to generate
s.type = 'QAM'; % subcarrier constellation type
s.SNRdB = SNRdb; % signal to noise ratio, in dB
s.Fsr = s.Fs; % simulated receiver sampling rate
s.Fcr = 11.805e9; % simulated receiver center frequency
s.Fc = 11.805e9; % OFDB Symbol center frequency
s.beta = 0; % CFO
s.gutter = 1; % boolean for gutters
s.data = randi([0 s.Midx-1],s.N,Nprovided); % Input data


% Generate frame, zero pad, and process
y = genOFDM(s);
y = [y; zeros((N+Ng)*ceil(length(y)/(N+Ng)) - length(y),1)]; % zero pad

p.N = N;
p.Ng = Ng;
p.Nsym = Nsym;
p.y = y;

[snr_dB, nitter] = ofdmSnrEstimator(p);



y = reshape(y,N+Ng,[]); % break up in symbols
y = y(Ng+1:end,:); % remove CP
Y = 1/sqrt(N).*fft(y);


gutter_check = 0;
if (abs(Y(1)) < epsilon && abs(Y(2)) < epsilon && ...
        abs(Y(end)) < epsilon && abs(Y(end-1)) < epsilon)
    gutter_check = 1;
end

% set gutter values to nan() for easier comparison
if (s.gutter)
    Y(1:2,:) = nan();
    Y(end-1:end,:) = nan();
end

% Expected decoded input data
if ( upper(string(s.type)) == "PSK")
    for ii=1:Nprovided
     seq = pskmod(s.data(:,ii),2^nextpow2(max(s.data(:,ii))));
     if (s.gutter)
     seq(1:2) = nan();
     seq(end-1:end) = nan();
     end
     data_in(:,ii) = getDiffEnc(seq);
    end
elseif( upper(string(s.type)) == "QAM")
    for ii=1:Nprovided
     seq = qammod(s.data(:,ii),2^nextpow2(max(s.data(:,ii))),'UnitAveragePower',true);
     if (s.gutter)
     seq(1:2) = nan();
     seq(end-1:end) = nan();
     end
     data_in(:,ii) = getDiffEnc(seq);
    end
end
if (Nsym - Nprovided > 0)
    Nreps = floor(Nsym/Nprovided);
    Nremain = mod(Nsym,Nprovided);
    data_in = repmat(data_in,[1,Nreps]);
    data_in = [data_in , data_in(:,1:Nremain)];
end

% decoded output data
data_out = getDiffEnc(Y);

% Compare
comp = data_out==data_in;
for ii = 1:Nsym
    match(ii) = sum(isnan(data_in(isnan(data_out(:,ii)),ii)));
end
match = sum(comp) + match;
if ( mean(match) > avgSCdecodethresh)
    str = 'pass';
else
    str = 'Failed';
    warning("Test failed, %.2f / %d symbols correctly decoded on average",mean(match),Nsym)
end
if (abs(snr_dB-SNRdb)>snrthresh_dB)
    str = 'Failed';
    warning("Test failed, SNR estimate off by %.2f dB",abs(snr_dB-SNRdb))
end
fprintf('4. OFDM Test Noise, with gutter, full BW, %d symbols : %16s \n',Nsym,str)

clear all;
% ----------- Test 5 
% No noise, no gutter, partial BW, 1 symbol

N = 1024; % # of s.c
Ng = 32; % c.p length
Midx = 4; % subcarrier constellation size
Nsym = 1; % # of symbols to simulate
Nprovided = 1; % # of provided data sequences to generate
SNRdb = nan();
Fsr = 100e6;

% Test input params
s.Fs = 240e6; % Symbol BW
s.N = N; % # of s.c.
s.Ng = Ng; % c.p length
s.Midx = Midx; % subcarrier constellation size
s.Nsym = Nsym; % # of consequtive frames to generate
s.type = 'QAM'; % subcarrier constellation type
s.SNRdB = SNRdb; % signal to noise ratio, in dB
s.Fsr = Fsr; % simulated receiver sampling rate
s.Fcr = 11.805e9; % simulated receiver center frequency
s.Fc = 11.805e9; % OFDB Symbol center frequency
s.beta = 0; % CFO
s.gutter = 1; % boolean for gutters
s.data = randi([0 s.Midx-1],s.N,Nprovided); % Input data


% Generate frame, zero pad, and process
y = genOFDM(s);
% Resample
if (s.Fsr < s.Fs)
%     tfcr = ((0:(length(y)-1))./s.Fsr)';
%     tfc = (0:1/s.Fs:tfcr(end))';
%     y = interp1(tfcr,y,tfc);
    tVec = (0:length(y)-1)'/s.Fsr;
    y = resample(y,tVec,s.Fs);
end

y = [y; zeros((N+Ng)*ceil(length(y)/(N+Ng)) - length(y),1)]; % zero pad
y = reshape(y,N+Ng,[]); % break up in symbols
y = y(Ng+1:end,:); % remove CP
Y = 1/sqrt(N).*fft(y);

% Expected decoded input data
if ( upper(string(s.type)) == "PSK")
    for ii=1:Nprovided
     seq = pskmod(s.data(:,ii),2^nextpow2(max(s.data(:,ii))));
     data_in(:,ii) = getDiffEnc(seq);
    end
elseif( upper(string(s.type)) == "QAM")
    for ii=1:Nprovided
     seq = qammod(s.data(:,ii),2^nextpow2(max(s.data(:,ii))),'UnitAveragePower',true);
     data_in(:,ii) = getDiffEnc(seq);
    end
end
if (Nsym - Nprovided > 0)
    Nreps = floor(Nsym/Nprovided);
    Nremain = mod(Nsym,Nprovided);
    data_in = repmat(data_in,[1,Nreps]);
    data_in = [data_in , data_in(:,1:Nremain)];
end

% decoded output data
data_out = getDiffEnc(Y);

% Compare
comp = data_out==data_in;
for ii = 1:Nsym
    match(ii) = sum(isnan(data_in(isnan(data_out(:,ii)),ii)));
end
match = sum(comp) + match;
str = 'pass';
if (floor(s.Fsr/s.Fs*N) > match)
    str = 'fail';
end
fprintf('5. OFDM Test no noise, no gutter, %.1f/%.1f MHz BW : %17s \n',s.Fsr/(1e6),s.Fs/(1e6),str)



%--------------------------------------%
%-----Starlink OFDM SIGNAL TEST--------%
%--------------------------------------%
clear all;
% ----------- Test 1 
% No noise, Full BW, 1 symbol

N = 1024; % # of s.c for starlink
Ng = 32; % c.p length for starlink
Midx = 4; % subcarrier constellation size
Nsym = 1; % # of symbols to simulate
Nprovided = 1; % # of provided data sequences to generate


% Test input params
s.Nsym = 1;
s.Midx = 4;
s.type = 'QAM';
s.SNRdB = nan();
s.Fsr = 240e6;
s.Fcr = 11.805e9;
F = 240e6/1024;
chIdx = round((s.Fcr/1e9 - 10.7 - F/2/1e9)/0.25 + 0.5);
Fcii = (10.7e9 + F/2 + 250e6*(chIdx - 0.5));
s.Fcr = Fcii; % Center at Starlink channel center
s.beta = 0; % CFO
s.data = randi([0 s.Midx-1],N,Nprovided); % Input data

% Generate frame, zero pad, and process
y = genStrlkOFDM(s);
y = [y; zeros((N+Ng)*ceil(length(y)/(N+Ng)) - length(y),1)]; % zero pad
y = reshape(y,N+Ng,[]); % break up in symbols
y = y(Ng+1:end,:); % remove CP
Y = 1/sqrt(N).*fft(y);

% set gutter values to nan() for easier comparison
Y(1:2,:) = nan();
Y(end-1:end,:) = nan();


% Expected decoded input data
if ( upper(string(s.type)) == "PSK")
    for ii=1:Nprovided
     seq = pskmod(s.data(:,ii),2^nextpow2(max(s.data(:,ii))));
     seq(1:2) = nan();
     seq(end-1:end) = nan();
     data_in(:,ii) = getDiffEnc(seq);
    end
elseif( upper(string(s.type)) == "QAM")
    for ii=1:Nprovided
     seq = qammod(s.data(:,ii),2^nextpow2(max(s.data(:,ii))),'UnitAveragePower',true);
     seq(1:2) = nan();
     seq(end-1:end) = nan();
     data_in(:,ii) = getDiffEnc(seq);
    end
end
if (Nsym - Nprovided > 0)
    Nreps = floor(Nsym/Nprovided);
    Nremain = mod(Nsym,Nprovided);
    data_in = repmat(data_in,[1,Nreps]);
    data_in = [data_in , data_in(:,1:Nremain)];
end

% decoded output data
data_out = getDiffEnc(Y);

% Compare
comp = data_out==data_in;
for ii = 1:Nsym
    match(ii) = sum(isnan(data_in(isnan(data_out(:,ii)),ii)));
end
match = sum(comp) + match;
if (Nsym == sum(match == N))
    str = 'pass';
else
    str = 'Failed';
    warning("Test failed, %d / %d symbols correctly decoded",match,N)
end
fprintf('5. Starlink symbol Test no noise, full BW, 1 symbol : %18s \n',str)

clear all;
% ----------- Test 2 
% No noise, Full BW, multiple symbol

N = 1024; % # of s.c for starlink
Ng = 32; % c.p length for starlink
Midx = 4; % subcarrier constellation size
Nsym = 300; % # of symbols to simulate
Nprovided = 1; % # of provided data sequences to generate


% Test input params
s.Nsym = Nsym;
s.Midx = 4;
s.type = 'QAM';
s.SNRdB = nan();
s.Fsr = 240e6;
s.Fcr = 11.805e9;
F = 240e6/1024;
chIdx = round((s.Fcr/1e9 - 10.7 - F/2/1e9)/0.25 + 0.5);
Fcii = (10.7e9 + F/2 + 250e6*(chIdx - 0.5));
s.Fcr = Fcii; % Center at Starlink channel center
s.beta = 0; % CFO
s.data = randi([0 s.Midx-1],N,Nprovided); % Input data

% Generate frame, zero pad, and process
y = genStrlkOFDM(s);
y = [y; zeros((N+Ng)*ceil(length(y)/(N+Ng)) - length(y),1)]; % zero pad
y = reshape(y,N+Ng,[]); % break up in symbols
y = y(Ng+1:end,:); % remove CP
Y = 1/sqrt(N).*fft(y);

% set gutter values to nan() for easier comparison
Y(1:2,:) = nan();
Y(end-1:end,:) = nan();


% Expected decoded input data
if ( upper(string(s.type)) == "PSK")
    for ii=1:Nprovided
     seq = pskmod(s.data(:,ii),2^nextpow2(max(s.data(:,ii))));
     seq(1:2) = nan();
     seq(end-1:end) = nan();
     data_in(:,ii) = getDiffEnc(seq);
    end
elseif( upper(string(s.type)) == "QAM")
    for ii=1:Nprovided
     seq = qammod(s.data(:,ii),2^nextpow2(max(s.data(:,ii))),'UnitAveragePower',true);
     seq(1:2) = nan();
     seq(end-1:end) = nan();
     data_in(:,ii) = getDiffEnc(seq);
    end
end
if (Nsym - Nprovided > 0)
    Nreps = floor(Nsym/Nprovided);
    Nremain = mod(Nsym,Nprovided);
    data_in = repmat(data_in,[1,Nreps]);
    data_in = [data_in , data_in(:,1:Nremain)];
end

% decoded output data
data_out = getDiffEnc(Y);

% Compare
comp = data_out==data_in;
for ii = 1:Nsym
    match(ii) = sum(isnan(data_in(isnan(data_out(:,ii)),ii)));
end
match = sum(comp) + match;
if (Nsym == sum(match == N))
    str = 'pass';
else
    str = 'Failed';
    warning("Test failed, %d / %d symbols correctly decoded",match,N)
end
fprintf('6. Starlink symbol Test no noise, full BW, %d symbols : %15s \n',Nsym,str)

clear all;
% ----------- Test 3 
% Noise, Full BW, multiple symbol

N = 1024; % # of s.c for starlink
Ng = 32; % c.p length for starlink
Midx = 4; % subcarrier constellation size
Nsym = 300; % # of symbols to simulate
Nprovided = 1; % # of provided data sequences to generate
SNRdb = 5;
avgSCdecodethresh = 800;
snrthresh_dB = 2;

% Test input params
s.Nsym = Nsym;
s.Midx = 4;
s.type = 'QAM';
s.SNRdB = SNRdb;
s.Fsr = 240e6;
s.Fcr = 11.805e9;
F = 240e6/1024;
chIdx = round((s.Fcr/1e9 - 10.7 - F/2/1e9)/0.25 + 0.5);
Fcii = (10.7e9 + F/2 + 250e6*(chIdx - 0.5));
s.Fcr = Fcii; % Center at Starlink channel center
s.beta = 0; % CFO
s.data = randi([0 s.Midx-1],N,Nprovided); % Input data

% Generate frame, zero pad, and process
y = genStrlkOFDM(s);
y = [y; zeros((N+Ng)*ceil(length(y)/(N+Ng)) - length(y),1)]; % zero pad

p.N = N;
p.Ng = Ng;
p.Nsym = Nsym;
p.y = y;

[snr_dB, nitter] = ofdmSnrEstimator(p);


y = reshape(y,N+Ng,[]); % break up in symbols
y = y(Ng+1:end,:); % remove CP
Y = 1/sqrt(N).*fft(y);

% set gutter values to nan() for easier comparison
Y(1:2,:) = nan();
Y(end-1:end,:) = nan();


% Expected decoded input data
if ( upper(string(s.type)) == "PSK")
    for ii=1:Nprovided
     seq = pskmod(s.data(:,ii),2^nextpow2(max(s.data(:,ii))));
     seq(1:2) = nan();
     seq(end-1:end) = nan();
     data_in(:,ii) = getDiffEnc(seq);
    end
elseif( upper(string(s.type)) == "QAM")
    for ii=1:Nprovided
     seq = qammod(s.data(:,ii),2^nextpow2(max(s.data(:,ii))),'UnitAveragePower',true);
     seq(1:2) = nan();
     seq(end-1:end) = nan();
     data_in(:,ii) = getDiffEnc(seq);
    end
end
if (Nsym - Nprovided > 0)
    Nreps = floor(Nsym/Nprovided);
    Nremain = mod(Nsym,Nprovided);
    data_in = repmat(data_in,[1,Nreps]);
    data_in = [data_in , data_in(:,1:Nremain)];
end

% decoded output data
data_out = getDiffEnc(Y);

% Compare
comp = data_out==data_in;
for ii = 1:Nsym
    match(ii) = sum(isnan(data_in(isnan(data_out(:,ii)),ii)));
end
match = sum(comp) + match;
if ( mean(match) > avgSCdecodethresh)
    str = 'pass';
else
    str = 'Failed';
    warning("Test failed, %.2f / %d symbols correctly decoded on average",mean(match),Nsym)
end
if (abs(snr_dB-SNRdb)>snrthresh_dB)
    str = 'Failed';
    warning("Test failed, SNR estimate off by %.2f dB",abs(snr_dB-SNRdb))
end
fprintf('7. Starlink symbol Test noise, full BW, %d symbols : %18s \n',Nsym,str)





