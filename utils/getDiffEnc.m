function s = getDiffEnc(x)
% getDifferentialEncoding returns the differentially decoded sequence x. 
%
% -- Input --
% x    N x K vector whose elements describe a differential 4QAM 
%      encoded sequence. If only Ns < N elements of the sequence are known, 
%      still provide x to be N long, with nan() in the place of unkown elements.
%      The 2nd dimension is to pararellize sequences
% 
% -- Output -- 
% 
% s   N x 1 vector of differentially decoded elements. If an element could
%     not be differentially decoded, the value it takes is nan().
% 

[N,K] = size(x);

aim1 = x(1:end-1,:);
ai = x(2:end,:);
di = ai.*conj(aim1);
thetai = atan2(imag(di),real(di));
s = round(thetai.*2/pi);
% si = 2 and si = -2 both correspond to pi, and so are equivalent.  By
% convention, we will set these to si = 2.
s(s == -2) = 2;
s(s == -1) = 3;
s = [nan(1,K); s];

end
