function data = fftAcqStrlk(s)
%  fftAcq performs an FFT acquisition based on the PSS and SSS.
%
% -- Input --
% A struct with the following fields:
% 
% * Fsr         Receiver sampling frequency in Hz
%
% * Fcr         Receiver center frequency in Hz
%  
% * fmax        Max Frequnecy for frequency search in Hz
%  
% * fstep       Frequnecy step for frequency search in Hz
%  
% * fstepFine   Frequnecy step for fine frequency search in Hz
%  
% * fmin        Min Frequnecy for frequency search in Hz
%  
% * y           Ns x 1 vector of baseband data, where Ns corresponds to
%               samples at least equal in length to a Starlink frame in time.
%   
% * Pfa        False Alarm probability used for the threshold to determine
%              if a signal is present for a PRN. The total acquisition false 
%              alarm probability. The Pfa for each search cell can be derived
%              from this.
%
%   -- Output -- 
%
%   A struct with the following fields:
%  
%   * grid      F x Nk acquisition grid where F is created from fmax, fmin,
%               and fstep, and Nk is the number of samples of a Starlink
%               frame sampled at 240 Msps.
%   
%   * tauvec    Nk x 1 vector of considered delay times of Starlink time of
%               arrival of the first frame with respect to the data input
%               beginning
%
%   * fdvec     F  x 1 vector of considered doppler frequencies
%
%   * fdfine    fine doppler frequency estimate from peak
%   
%   * tau       time offset estimate from peak
%
%   * Nc        Length of local replica
%
%   * sigma2_c  Variance of local replica (non-zero samples)
%
%   * pkVec     Local replica used for acquisition.



debugPltHT_enable = 0;
debugPltGrid_en = 0;
fontsize = 12;

%------------------------------------%
%----------- Params------------------%
%------------------------------------%
N = 1024; % # of sc
Ng = 32;  % CP length
Ns = N + Ng; % # samples per symbol
Ndsym = 300; % # of data-carying symbols
Fs = 240e6; % BW
Nk = 1/750*Fs; % # samples per frame
PSS = genPss();
SSS = genSss();
c = [PSS;SSS]; % known frame
pkVec = c;
c = [c; zeros(Nk-length(c),1)];
%------------------------------------%
%-------- Resample to full BW--------%
%------------------------------------%
y = s.y;
if (s.Fsr ~= Fs)
     tVec = [0:length(y)-1]'/(s.Fsr); 
     y = resample(y,tVec,Fs);
end
buffer = 100;
if (length(y) >= Nk - buffer)
    y = [y;zeros(Nk-length(y),1)];
elseif (length(y) - Nk < buffer)
    str = sprintf("Error! The data should be long enough to contain at least 1 frame.");
    error(str)
end
%------------------------------------%
%-------- Remove receiver bias-------%
%------------------------------------%
if (getClosestFch(s.Fcr) - s.Fcr ~= 0)
    tVec = [0:length(y)-1]'/(Fs); 
    Fshift =  getClosestFch(s.Fcr) - s.Fcr;
    y = y.*exp(-1j*2*pi*Fshift*tVec);
end

%------------------------------------%
%-----------Derived Params-----------%
%------------------------------------%
fdvec = (s.fmin : s.fstep : s.fmax)';
tauvec = (0:(Nk-1))'/Fs;
known = c(c ~= 0);
Nc = length(known);
sigma2_c = var(known);
Ngrid = length(fdvec).*length(tauvec); % Number of cells in the grid
C = fft(c);
Nac = floor(length(y)/Nk);
y = y(1:(Nac*Nk)); % Multiple of frame length
% Initializations
values = zeros(4,1); % will hold [fdfine,tau,SNRdb,detection]
grids = zeros(length(fdvec),Nk);
Sk = zeros(length(fdvec),length(c), Nac);
%------------------------------------%
%----------Generate Acq. Grid--------%
%------------------------------------%
for jj = 1:Nac
    for ii = 1:length(fdvec)
        pre_idcs = ((jj-1)*Nk+1):(jj*Nk);
        tau_pre = tauvec(1:length(pre_idcs));
        x = y(pre_idcs);
        beta = -fdvec(ii)/getClosestFch(s.Fcr);
        FsD = (1 + beta)*Fs;
        [x,tau_post] = resample(x,tau_pre,FsD);
        
        th_hat = 2*pi*(fdvec(ii))*(tau_post);
        
        x_tilde = x.*exp(-j*th_hat);
        x_tilde = [x_tilde; zeros(Nk-length(x_tilde),1)];
        x_tilde = x_tilde(1:Nk);

        X_tilde = fft(x_tilde);
        Zr = X_tilde.*conj(C);
        sk = ifft(Zr)';
        Sk(ii,:,jj) = sk(1:length(c));
    end
end 
grids(:,:) = sum(abs(Sk),3); % Save grid

%------------------------------------%
%----------Process Acqiosition-------%
%------------------------------------%

% Find Peak, tau, and doppler
[tau_prob, fidx] = max(grids); % find max of frq values
[mx_val, tau_idx] = max(tau_prob); % find max of tau values
tau = tauvec(tau_idx);
fdcoarse = fdvec(fidx(tau_idx));

% Find finer doppler
fdfine = fdcoarse;
expSk2_fine = mx_val./Nac;
if ((s.fstep ~= 1)) % Finer look if step is not 1
    fdvec_fine = [(fdcoarse-floor(s.fstep/2)):s.fstepFine:(fdcoarse+floor(s.fstep/2))]';
    % Each column is a different phase shift
    thmat = 2*pi*(tauvec-tauvec(Nk))*(fdvec_fine)';
    z_fine = zeros(length(fdvec_fine),1);
    ctau = circshift(c,tau_idx-1);
    for jj=1:Nac
        iidum = ((jj-1)*Nk+1):(jj*Nk);
        x1_shift = y(iidum).*exp(-1j*thmat(1:length(iidum),:));
        z_fine(:) = z_fine(:) + abs(x1_shift'*ctau);     
    end
    [mx,fd_fine_idx] = max(z_fine);
    expSk2_fine = mx./Nac;
    fdfine = fdvec_fine(fd_fine_idx);
end

%------------------------------------%
%---------   Output -----------------%
%------------------------------------%

data.grid = grids;
data.tauvec = tauvec;
data.fdvec = fdvec;
data.fdfine = fdfine;
data.tau = tau;
data.Nc = Nc;
data.sigma2_c = sigma2_c; 
data.pkVec = pkVec;

data.Sk = Sk;

if (debugPltGrid_en)
    [X,Y] = meshgrid(tauvec,fdvec);
    figure()
    surf(X,Y,grids(:,:),'linestyle','none','FaceColor','interp','EdgeColor','none','FaceLighting','gouraud')
end
end
