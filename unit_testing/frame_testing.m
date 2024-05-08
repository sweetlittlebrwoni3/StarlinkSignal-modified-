clc;
clear all;
addpath('../Utils')

%--------------------------------------%
%--------Starlink FRAME TEST-----------%
%--------------------------------------%

% ----------- Test 1

N = 1024; % # of s.c
Ng = 32; % c.p length
Midx = 4; % subcarrier constellation size
Nsym = 300; % # of symbols to simulate
Nprovided = Nsym; % # of provided data sequences to generate

% Test input params
s.SNRdB = nan(); % Signal to noise ratio, in dB
s.Fsr = 240e6; % Receiver sample rate, in Hz. 
s.Fcr = 12075117187.5; % Receiver center frequency, in Hz
s.beta = 0; % Doppler parameter
s.type = 'QAM'; % subcarrier constellation type
s.data = randi([0 Midx-1],N,Nprovided); % Input data

% Generate frame, zero pad, and process
y = genStrlkFrame(s);
y = [y; zeros((N+Ng)*ceil(length(y)/(N+Ng)) - length(y),1)]; % zero pad
y = reshape(y,N+Ng,[]); % break up in symbols
y = y(Ng+1:end,3:302); % remove CP
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
    warning("Test failed, %d / %d symbols correctly decoded",sum(match == N),Nsym)
end
fprintf('1. Starlink Frame Test no noise : %30s \n',str)

clear all;
% ----------- Test 2

N = 1024; % # of s.c
Ng = 32; % c.p length
Midx = 4; % subcarrier constellation size
Nsym = 300; % # of data symbols in starlink frame
Nprovided = 1; % # of provided data sequences to generate
SNRdb = 5;
avgSCdecodethresh = 0.8;
snrthresh_dB = 2;

% Test input params
s.SNRdB = SNRdb; % Signal to noise ratio, in dB
s.Fsr = 240e6; % Receiver sample rate, in Hz. 
s.Fcr = 12075117187.5; % Receiver center frequency, in Hz
s.beta = 0; % Doppler parameter
s.type = 'QAM'; % subcarrier constellation type
s.data = randi([0 Midx-1],N,Nprovided); % Input data

% Generate frame, zero pad, and process
y = genStrlkFrame(s);
y = [y; zeros((N+Ng)*ceil(length(y)/(N+Ng)) - length(y),1)]; % zero pad

y = reshape(y,N+Ng,[]); % break up in symbols

p.N = N;
p.Ng = Ng;
p.Nsym = Nsym;
datasyms = y(:,3:302);
p.y = datasyms;
[snr_dB, nitter] = ofdmSnrEstimator(p);

y = y(Ng+1:end,3:302); % remove CP & PSS and SSS
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
if ( mean(match)/N > avgSCdecodethresh)
    str = 'pass';
else
    str = 'Failed';
    warning("Test failed, %.2f / %d subcarriers correctly decoded on average",mean(match),N)
end
if (abs(snr_dB-SNRdb)>snrthresh_dB)
    str = 'Failed';
    warning("Test failed, SNR estimate off by %.2f dB",abs(snr_dB-SNRdb))
end
fprintf('2. Starlink Frame Test noise, full BW : %24s \n',str)


clear all;
% ----------- Test 3 
% No noise, , partial BW, 1 frame

N = 1024; % # of s.c
Ng = 32; % c.p length
Fs = 240e6; % Starlink signal BW
Midx = 4; % subcarrier constellation size
Nsym = 300; % # of symbols to simulate
Nprovided = Nsym; % # of provided data sequences to generate
Fsr = 120e6;

% Test input params
s.SNRdB = nan(); % Signal to noise ratio, in dB
s.Fsr = Fsr; % Receiver sample rate, in Hz. 
s.Fcr = 12075117187.5; % Receiver center frequency, in Hz
s.beta = 0; % Doppler parameter
s.type = 'QAM'; % subcarrier constellation type
s.data = randi([0 Midx-1],N,Nprovided); % Input data

% Generate frame, resample, and process
y = genStrlkFrame(s);
% Resample
if (s.Fsr < Fs)
    tVec = (0:length(y)-1)'/s.Fsr;
    y = resample(y,tVec,Fs);
end
y = [y; zeros((N+Ng)*ceil(length(y)/(N+Ng)) - length(y),1)]; % zero pad
y = reshape(y,N+Ng,[]); % break up in symbols
y = y(Ng+1:end,:); % remove CP
Y = 1/sqrt(N).*fft(y);

% Disregard PSS and SSS for decoding
Y = Y(:,3:302);

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
if (floor(s.Fsr/Fs*N) > match)
    str = 'fail';
end
fprintf('3. Starlink Frame no noise, %.1f/%.1f MHz BW : %15s \n',s.Fsr/(1e6),Fs/(1e6),str)

clear all;
% ----------- Test 4 
% noise, , full BW, 1 frame, doppler
smallgrid = 1;
N = 1024; % # of s.c
Ng = 32; % c.p length
Fs = 240e6; % Starlink signal BW
Midx = 4; % subcarrier constellation size
Nsym = 300; % # of symbols to simulate
Nprovided = Nsym; % # of provided data sequences to generate
Fd_target = 200e3;
Fcr = getClosestFch(12075117187);
beta = -Fd_target./getClosestFch(Fcr);
offset = 200e3;
SNRdb = 15;
if (smallgrid)
    fmax = Fd_target+offset;
    fstep = 5e3; 
    fmin = Fd_target-offset;
else
    fmax = 400e3;
    fstep = 5e3; 
    fmin = -400e3;
end
fstepFine = 10;
Pfa = 0.01;
threshold_dopp = 500;
threshold_tau = 1/Fs;
tau_start = 0;
avgSCdecodethresh = 0.95;
snrthresh_dB = 2;

% Test input params
s.SNRdB = SNRdb; % Signal to noise ratio, in dB
s.Fsr = Fs; % Receiver sample rate, in Hz. 
s.Fcr = Fcr; % Receiver center frequency, in Hz
s.beta = beta; % Doppler parameter
s.type = 'QAM'; % subcarrier constellation type
s.data = randi([0 Midx-1],N,Nprovided); % Input data


% Generate frame, resample, and process
y = genStrlkFrame(s);
% Extract doppler estimate
a.Fsr = s.Fsr;
a.Fcr = s.Fcr;
a.fmax = fmax;
a.fstep = fstep;
a.fstepFine = fstepFine;
a.fmin = fmin;
a.y = y;
a.Pfa = Pfa;
out = fftAcqStrlk(a);

% Resample
if (s.Fsr < Fs)
    tVec = (0:length(y)-1)'/s.Fsr;
    y = resample(y,tVec,Fs);
end
y = [y; zeros((N+Ng)*ceil(length(y)/(N+Ng)) - length(y),1)]; % zero pad

str = 'pass';
if (out.values(4) ~= 1)
    str = 'fail';
    warning("Frame incorectly not detected")
elseif (abs(out.values(1)-Fd_target) > threshold_dopp)
    str = 'fail';
    warning("Doppler off by %.1f",abs(out.values(1)-Fd_target))
elseif (abs(out.values(2)-tau_start) > threshold_tau)
    str = 'fail';
    warning("Tau off by %.1f",abs(out.values(2)-tau_start))
elseif (~isnan(s.SNRdB) && abs(out.values(3)-SNRdb) > snrthresh_dB)
    str = 'fail';
    warning("SNR estimate off by %.1f",abs(out.values(3)-SNRdb))
end

% Doppler correct
tVec = [0:length(y)-1]'/(Fs); 
Fshift =  getClosestFch(s.Fcr) - s.Fcr + out.values(1);
y = y.*exp(-1j*2*pi*Fshift*tVec);
y = reshape(y,N+Ng,[]); % break up in symbols
y = y(Ng+1:end,:); % remove CP
Y = 1/sqrt(N).*fft(y);

% Disregard PSS and SSS for decoding
Y = Y(:,3:302);

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
if ( mean(match)/N > avgSCdecodethresh)
    str = 'pass';
else
    str = 'Failed';
    warning("Test failed, %.2f / %d subcarriers correctly decoded on average",mean(match),N)
end
fprintf('4. Starlink Frame noise, doppler, %.1f/%.1f MHz BW : %9s \n',s.Fsr/(1e6),Fs/(1e6),str)


clear all;
% ----------- Test 5 
% No noise, , full BW, 1 frame, doppler
smallgrid = 1;
N = 1024; % # of s.c
Ng = 32; % c.p length
Fs = 240e6; % Starlink signal BW
Fsr = 120e6;
Midx = 4; % subcarrier constellation size
Nsym = 300; % # of symbols to simulate
Nprovided = Nsym; % # of provided data sequences to generate
Fd_target = 200e3;
Fcr = getClosestFch(12075117187);
beta = -Fd_target./getClosestFch(Fcr);
offset = 200e3;
SNRdb = 15;
if (smallgrid)
    fmax = Fd_target+offset;
    fstep = 5e3; 
    fmin = Fd_target-offset;
else
    fmax = 400e3;
    fstep = 5e3; 
    fmin = -400e3;
end
fstepFine = 10;
Pfa = 0.01;
threshold_dopp = 500;
threshold_tau = 1/Fs;
tau_start = 0;
avgSCdecodethresh = 0.95;
snrthresh_dB = 2.5;

% Test input params
s.SNRdB = SNRdb; % Signal to noise ratio, in dB
s.Fsr = Fsr; % Receiver sample rate, in Hz. 
s.Fcr = Fcr; % Receiver center frequency, in Hz
s.beta = beta; % Doppler parameter
s.type = 'QAM'; % subcarrier constellation type
s.data = randi([0 Midx-1],N,Nprovided); % Input data


% Generate frame, resample, and process
y = genStrlkFrame(s);
% Extract doppler estimate
a.Fsr = s.Fsr;
a.Fcr = s.Fcr;
a.fmax = fmax;
a.fstep = fstep;
a.fstepFine = fstepFine;
a.fmin = fmin;
a.y = y;
a.Pfa = Pfa;
out = fftAcqStrlk(a);

% Resample
if (s.Fsr < Fs)
    tVec = (0:length(y)-1)'/s.Fsr;
    y = resample(y,tVec,Fs);
end
y = [y; zeros((N+Ng)*ceil(length(y)/(N+Ng)) - length(y),1)]; % zero pad


str = 'pass';
if (out.values(4) ~= 1)
    str = 'fail';
    warning("Frame incorectly not detected")
elseif (abs(out.values(1)-Fd_target) > threshold_dopp)
    str = 'fail';
    warning("Doppler off by %.1f",abs(out.values(1)-Fd_target))
elseif (abs(out.values(2)-tau_start) > threshold_tau)
    str = 'fail';
    warning("Tau off by %.1f",abs(out.values(2)-tau_start))
elseif (~isnan(s.SNRdB) && abs(out.values(3)-SNRdb) > snrthresh_dB)
    str = 'fail';
    warning("SNR estimate off by %.1f",abs(out.values(3)-SNRdb))
end

% Doppler correct
tVec = [0:length(y)-1]'/(Fs); 
Fshift =  getClosestFch(s.Fcr) - s.Fcr + out.values(1);
y = y.*exp(-1j*2*pi*Fshift*tVec);
y = reshape(y,N+Ng,[]); % break up in symbols
y = y(Ng+1:end,:); % remove CP
Y = 1/sqrt(N).*fft(y);

% Disregard PSS and SSS for decoding
Y = Y(:,3:302);

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

% Only look at valid indices
data_out = data_out(getValidSC(s.Fcr,s.Fsr),:);
data_in = data_in(getValidSC(s.Fcr,s.Fsr),:);
Nvalid = length(data_in);
% Compare
comp = data_out==data_in;
for ii = 1:Nsym
    match(ii) = sum(isnan(data_in(isnan(data_out(:,ii)),ii)));
end
match = sum(comp) + match;
if ( mean(match)/Nvalid > avgSCdecodethresh)
    str = 'pass';
else
    str = 'Failed';
    warning("Test failed, %.2f / %d subcarriers correctly decoded on average",mean(match),Nvalid)
end
fprintf('5. Starlink Frame noise, doppler, %.1f/%.1f MHz BW : %9s \n',s.Fsr/(1e6),Fs/(1e6),str)


