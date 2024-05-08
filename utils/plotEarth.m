function plotEarth(fig, rx_lla)
% Generates a 3D figure of Earth given a figure object
% 
%  -- Input --
%   fig     Figure object (i.e fig = figure())
%
%   rx_lla  (Optional) [lattitude, longitude, altitude]
%
addpath(fullfile(matlabroot,'examples/map/data'));

figure(fig)
grs80 = referenceEllipsoid("grs80","km");
ax = axesm("globe","Geoid",grs80);
ax.Position = [0 0 1 1];
axis equal off
view(3)
load topo60c
geoshow(topo60c,topo60cR,"DisplayType","texturemap")
demcmap(topo60c)
load coastlines
geoshow(coastlat,coastlon,"Color","black")
rivers = readgeotable("worldrivers.shp");
geoshow(rivers,"Color","blue")
if exist('rx_lla','var')
    rx_ECEF =  lla2ecef(rx_lla)'./1e3; % in km
    geoshow(rx_lla(1), rx_lla(2), 'DisplayType', 'Point', 'Marker', '+', 'Color', 'red');
    h = 500; % STarlink at ~ 500Km
n = 50; % number of points
r = h;
theta = linspace(0,2*pi,n); % angle around the z-axis
x = linspace(0,h,n); % height along the x-axis
[Theta,X] = meshgrid(theta,x); % create a grid of theta and z values
R = r*X/h; % calculate the radius at each height
Z = R.*cos(Theta); % convert to cartesian coordinates
Y = R.*sin(Theta);
% Define the rotation angles in degrees
pitch = 0; 
roll = rx_lla(1); % +lattitude
yaw = -rx_lla(2); % -longitude
% Convert the angles to radians
pitch = pitch*pi/180;
roll = roll*pi/180;
yaw = yaw*pi/180;
% Define the rotation matrices
Rx = [1 0 0; 0 cos(pitch) -sin(pitch); 0 sin(pitch) cos(pitch)]; % rotation matrix around x-axis
Ry = [cos(roll) 0 sin(roll); 0 1 0; -sin(roll) 0 cos(roll)]; % rotation matrix around y-axis
Rz = [cos(yaw) -sin(yaw) 0; sin(yaw) cos(yaw) 0; 0 0 1]; % rotation matrix around z-axis
% Apply the rotations to the coordinates
XYZ = [X(:) Y(:) Z(:)]; % reshape the coordinates into a matrix of size (n^2,3)
XYZ = XYZ*Rx*Ry*Rz; % multiply by the rotation matrices
X = reshape(XYZ(:,1),n,n); % reshape back into a matrix of size (n,n)
Y = reshape(XYZ(:,2),n,n);
Z = reshape(XYZ(:,3),n,n);
% Offset to RX location
X = X+rx_ECEF(1); % reshape back into a matrix of size (n,n)
Y = Y+rx_ECEF(2);
Z = Z+rx_ECEF(3);
% Plot the cone
surf(X,Y,Z) % create a surface plot
end

end