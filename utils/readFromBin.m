function [y] = readFromBin(filepath, Fsr, tDur, tSeek)
% readFromBin returns the complex IQ samples from a bin file. 
%
% -- Input --
% filepath    A string of the path to the file.
%
% Fsr         The sampling rate of the data.
%
% tDur        The duration to read from the file in seconds.
%
% tSeek       The starting location in seconds to read the data from in the
%             file.

% -- Output -- 
% 
% y   (Tdur*Fsr) x 1 vector of complex samples
% 

fid = fopen(filepath, 'r', 'n'); % Read only, native byte ordering
Ns = floor(tDur*Fsr);
% 4 bytes per complex sample
seekOffset = floor(tSeek*Fsr)*4;
status = fseek(fid,seekOffset,-1);
if(status == -1)
    error('tSeek beyond file limit');
end

x = fread(fid, [2,Ns], 'int16')';
fclose(fid);
if(length(x(:,1)) < Ns)
   error('Insufficient data');
end

y = complex(x(:,1),x(:,2));

end