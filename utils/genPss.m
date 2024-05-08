function [pss] = genPss()
% genPss generates the Starlink downlink primary synchronization sequence (PSS). 
%
% -- Output -- 
%
%  pss    (N + Ng) x 1 PSS. For starlink N = 1024, Ng = 32. Returned
%         sampled at 240e6 Hz
%
% -- Author -- 
%  Dr. Todd Humphreys              
%   

% Duration of OFDM symbol guard interval (cyclic prefix), expressed in
% intervals of 1/Fs.  
Ng = 32;
% Number of sub-segments in the primary synchronization sequence (PSS)
NpssSeg = 8;

% Fibonacci LFSR
ciVec = [3 7]';
a0Vec = [0 0 1 1 0 1 0]';
n = 7;
m = 2^n - 1;
lfsrSeq = zeros(m,1);
for idx=1:m
    buffer = a0Vec(ciVec);
    val = rem(sum(buffer),2); % Cascaded XOR is the same as the parity
    a0Vec = [val; a0Vec(1:end-1)];
    lfsrSeq(idx) = val;
end
lfsrSeqMod = flipud(lfsrSeq);
lfsrSeqMod = [0; lfsrSeqMod];
seq = 2*lfsrSeqMod-1;
pkSegVec = [exp(-1j*pi/4 - 1j*0.5*pi*cumsum(seq))];
pkVec = [-pkSegVec; repmat(pkSegVec,NpssSeg - 1,1)];
pkCP = pkVec(end-Ng+1:end);
pss = [-pkCP; pkVec];

end