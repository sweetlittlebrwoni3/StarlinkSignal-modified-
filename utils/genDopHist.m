function [fdop,epoch_datetimes, epoch_state_interp] = genDopHist(s)
%  genDopHist generates a doppler history interpolated from the provided
%             parameters
%
% -- Input --
% A struct with the following fields:
% 
% * epoch_datetime   Nx1 vector of datetimes
%
% * epoch_state      Nx6 corresponding ECEF state [x,y,z,vx,vy,vz] in km
%
% * start_datetime   datime corresponding to beginning of desired simulated
%                    doppler history
%
% * tdur             duration of simulated doppler history
%
% * Fc               center frequency of signal
%
% * Fs               sampling frequency of desired simulated doppler history
%
% * rx_state         1x6 corresponding ECEF state of receiver
%                    [x,y,z,vx,vy,vz] in km
%  
%   -- Output -- 
%
% * fdop             A 

Nmax = 32000000;
c = physconst('Lightspeed')./(1e3); % in km/s
lamda = c./s.Fc; % in Km
Nsamples = s.Fs.*s.tdur;
% if (Nsamples > Nmax)
%     Nitter = ceil(Nsamples/Nmax);
%     fdop = tall(zeros(Nsamples,1));
%     epoch_datetimes = tall(zeros(Nsamples,1));
%     epoch_state_interp = tall(zeros(6,Nsamples));
%     
%     for ii=1:Nitter
%         Nsrt_i = (ii-1)*Nmax+1;
%         Nend_i = min(Nsrt_i + Nmax, Nsamples);
%         Nsamples_i = Nend_i-Nsrt_i;
%         strDateTime_i = s.start_datetime + seconds((Nsrt_i-1)./s.Fs);
% 
%         range = seconds([0:(Nsamples_i-1)]'./s.Fs);
% 
%         epoch_datetimes_i = strDateTime_i + range;
%         epoch_state_interp_i = zeros(6,length(epoch_datetimes_i));
%         for jj = 1:6
%             epoch_state_interp_i(ii,:) = interp1(s.epoch_datetime,s.epoch_state(ii,:),epoch_datetimes_i);
%         end
%         % Calculate the relative velocity between the satellite and the receiver
%         d = vecnorm(s.rx_state(1:3)-epoch_state_interp_i(1:3,:));
%         rG = (s.rx_state(1:3)-epoch_state_interp_i(1:3,:))./d; % Unit vector from SV to RX
%         
%         v = s.rx_state(4:6)-epoch_state_interp_i(4:6,:); % Relative velocity from SV to Rx
%         
%         % Calculate the doppler shift using the formula df = fapp - fc
%         fdop_i = - 1 / lamda * dot(rG,v)';
%         
%          fdop(Nsrt_i:Nend_i) = fdop_i;
%          epoch_datetimes(Nsrt_i:Nend_i) = epoch_datetimes_i;
%          epoch_state_interp(Nsrt_i:Nend_i) = epoch_state_interp_i;
%     end
% else 
    range = seconds([0:(Nsamples-1)]'./s.Fs);
        
    epoch_datetimes = s.start_datetime + range;
    epoch_state_interp = zeros(6,length(epoch_datetimes));
    for ii = 1:6
        epoch_state_interp(ii,:) = interp1(s.epoch_datetime,s.epoch_state(ii,:),epoch_datetimes);
    end
    
    % Calculate the relative velocity between the satellite and the receiver
    d = vecnorm(s.rx_state(1:3)-epoch_state_interp(1:3,:));
    rG = (s.rx_state(1:3)-epoch_state_interp(1:3,:))./d; % Unit vector from SV to RX
    
    v = s.rx_state(4:6)-epoch_state_interp(4:6,:); % Relative velocity from SV to Rx
    
    % Calculate the doppler shift using the formula df = fapp - fc
    fdop = - 1 / lamda * dot(rG,v)';
% end
end


