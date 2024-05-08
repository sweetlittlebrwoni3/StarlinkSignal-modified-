function [snr_dB, sigma_2, nitter] = ofdmSnrEstimator(s)
% ofdmSnrEstimator  estimates the SNR from a stream of OFDM symbols.
%                   The implementation is based off "The Two-Step 
%                   Non-Data-Aided SNR Estimation in the Low SNR
%                   Region of OFDM Signals" by D. WANG and W. XU.
%                   In regards to a starlink
% -- Input --
% A struct with the following fields:
% 
% * N               Number of subcarriers
% 
% * Ng              Length of cyclic prefix (CP)
% 
% * Nsym            Number of consequtive symbols in y
% 
% * y               Nsym*(N + Ng) x 1 vector of samples of the OFDM symbols
%                   in time
% 
% * threshold       (Optional) threshold dictating the ML itereration decision
%
%   -- Output -- 
%
%  snr_dB       Estimate of SNR in dB
%
%  sigma_2      Noise variance estimate
%
%  nitter       Number of itterations in ML to meet threshold

if (~isfield(s,'threshold'))
    s.threshold = 1e-6;
end
y = reshape(s.y,s.N+s.Ng,[]); % break up in symbols


u = 1:s.Ng;
Ju = zeros(length(u),1);
ksi = zeros(length(u),1);
Lhat = zeros(length(u),1);
for uu=length(u):-1:1
    Ju(uu) = 1/(2*s.Nsym*(s.Ng-uu+1)).*sum(sum(abs(y(uu:s.Ng,:)-y((s.N+uu):(s.N+s.Ng),:)).^2));
    if (uu == length(u))
        ksi(uu) = Ju(uu);
    else
        ksi(uu) = Ju(uu) - (1-1/(s.Ng-uu+1)).*Ju(uu); 
    end
    % given L = uu
    sigma_2 = Ju(uu);
    L_hat(uu) = prod(normpdf(ksi(uu:end),sigma_2/(s.Ng-uu+1),sigma_2^2/(s.Nsym.*(s.Ng-uu+1).^2))).^(1/(s.Ng-uu+1));
end
[L_hat, idx] = min(L_hat);
sigma_2_hat = Ju(idx);
s_hat = 1/(s.Nsym*s.N).*sum(sum(real(y(1:(s.Ng+Lhat),:).*conj(y(s.N+1:(s.Ng+Lhat+s.N),:)))));
Ey_hat = 1/(s.Nsym*(s.N-s.Ng)).*sum(sum(abs(y).^2));

theta_c_hat = [s_hat; sigma_2_hat; Ey_hat];
nitter = 0;
while (1)
    varShat = (s.Ng+L_hat).*(Ey_hat^2 + s_hat^2)./(2*s.Nsym*s.Ng^2);
    varsigma2hat = sigma_2_hat.^2./(s.Nsym.*(s.Ng-L_hat));
    varEyhat = Ey_hat^2./(s.Nsym.*(s.N-L_hat));
    covShatsigma2hat = -sigma_2_hat.^2/(2*s.Nsym*s.Ng);
    C = [varShat, covShatsigma2hat, 0;
        covShatsigma2hat, varsigma2hat, 0;
        0 , 0, varEyhat];
    H = [1 0; 0 1; 1 1];
    thetar = inv(H'*inv(C)*H)*H'*inv(C)*theta_c_hat;
    e = abs(thetar-theta_c_hat(1:2)).^2;
    theta_c_hat(1:2) = thetar;
    nitter = nitter + 1;
    if (sum(e < s.threshold) == 2)
        break;
    end

end
sigma_2 = thetar(2);
SNR_hat = thetar(1)./thetar(2);
snr_dB = 10*log10(SNR_hat);


end