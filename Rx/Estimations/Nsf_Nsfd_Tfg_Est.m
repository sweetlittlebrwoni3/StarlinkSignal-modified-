function [Nsf , Nsfd , Tfg] = Nsf_Nsfd_Tfg_Est(s)
% Nsf_Nsfd_Tfg_Est          Calculates the Nsf , Nsfd & Tfg parameters
% 
% -- Input --
% * N                       Number of subcarriers
%
% * Ng                      Cyclic prefix length
%
% * Fs                      Sampling frequency in the receiver
%
% * Tf                      Frame length
%
%   -- Output -- 
% * Nsf                     Number of non-zero symbols per frame
%
% * Nsfd                    Number of data symbols in a frame
%
% * Tfg                     Frame guard interval

N = s.N;
Ng = s.Ng;
Fs = s.Fs;
Tf = s.Tf;

Tsym = (N + Ng)/Fs;
Nsf = floor(Tf/Tsym) - 1;

Nsfd = Nsf - 4;

Tfg = Tf - Nsf*Tsym;

% s.Nsf = Nsf;
% s.Nsfd = Nsfd;
% s.Tfg = Tfg;

end