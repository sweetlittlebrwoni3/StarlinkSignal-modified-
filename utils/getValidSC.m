function idxs = getValidSC(Fcr, Fsr)
% getValidSc returns the valid subcarrier indices based on the provided 
%            center and sampling frequencies
%
%           |<---          Fs = 240 Mhz           --->|
%
%             $<- Fsr->$ 
%           -------------------------------------------
%           | $        $                              |
%           | $        $                              |
%  <--------------|----------------|-----------------------------> f (Hz)
%                Fcr              Fcii  
%                                  | 
%  <--------513-514-...-------1024-1-2-3-...----------512---------> s.c. idxs
%
% --- Input ---
%
% Fcr  receiver center frequency
% 
% Fsr  receiver sampling frequency
%

Fs = 240e6;
F = Fs/1024; % Starlink OFDM subcarrier BW

% Get closest channel center
chIdx = round((Fcr/1e9 - 10.7 - F/2/1e9)/0.25 + 0.5);
Fcii = (10.7e9 + F/2 + 250e6*(chIdx - 0.5));
% Find s.c. channel centers
Fsc = (Fcii-round(Fs/2-F/2)):F:(Fcii+round(Fs/2-F/2));
% Find capture frequency range
Fcr_start = Fcr - round(Fsr/2); 
Fcr_end = Fcr + round(Fsr/2);
% Find s.c. index
Fsc_start_idx = (Fsc - F/2 - Fcr_start);
Fsc_start_idx = find(Fsc_start_idx > 0 ,1);
if (isempty(Fsc_start_idx))
    error("Frequency range does not match as beginning of band is outside a starlink channel.")
end
Fsc_end_idx = (Fsc + F/2 - Fcr_end);
Fsc_end_idx = min(find(Fsc_end_idx > -1 ,1),1024);
if (isempty(Fsc_end_idx))
    Fsc_end_idx = 1024;
    warning("Frequency range goes beyond the channel. Truncating output to edge of channel.")
end
% Map 1-1024 to [513:1024,1:512]
idxs = Fsc_start_idx:Fsc_end_idx;
idxs = mod(idxs+511,1024)+1;

end
