function Fci = Fci_Est(s)
% Fci_Est                         Estimates the center frequency of ith
%                                 channel
% 
% -- Input --
% * data                          Input signal
%
% * N                             Number of subcarriers
%
% * Ng                            Cyclic prefix length
%
% * Fs                            Sampling frequency in the receiver
%
% * bs                            Bits per symbol
%
% * FcEst                         A priori estimate of Fc and
%                                 the exact center of the band captured to 
%                                 produce the input
%
%   -- Output -- 
% * Fci                           The value for maximum comparison
%







result = round((FcEst)/(1 + betam0Est - betam0));

Fci = result;

end