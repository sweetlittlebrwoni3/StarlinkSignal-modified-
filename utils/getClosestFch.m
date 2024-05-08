function Fcii = getClosestFch(Fc)
%  getClosestFch returns the center frequency of the closses Starlink
%                OFDM channel
% --- Input ---
%
% Fc    Center frequency of receiver in Hz
%
F = 240e6/1024;
chIdx = round((Fc/1e9 - 10.7 - F/2/1e9)/0.25 + 0.5);
Fcii = (10.7e9 + F/2 + 250e6*(chIdx - 0.5));

end