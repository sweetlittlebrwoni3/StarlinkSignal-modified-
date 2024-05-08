function [fdop,epoch_datetimes, epoch_state_interp] = genPLHist(s)
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
% * exp              Path loss exponent
%  
%   -- Output -- 
%
% * fdop             A 

    Nmax = 32000000;
    c = physconst('Lightspeed')./(1e3); % in km/s
    lamda = c./s.Fc; % in Km
    Nsamples = s.Fs.*s.tdur;
    range = seconds([0:(Nsamples-1)]'./s.Fs);
        
    epoch_datetimes = s.start_datetime + range;
    epoch_state_interp = zeros(6,length(epoch_datetimes));
    for ii = 1:6
        epoch_state_interp(ii,:) = interp1(s.epoch_datetime,s.epoch_state(ii,:),epoch_datetimes);
    end
    
    % Calculate the distance velocity between the satellite and the receiver
    d = vecnorm(s.rx_state(1:3)-epoch_state_interp(1:3,:));
    


    
end


