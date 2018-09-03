function Rg = sensorSurfaceGrid(Rs,Ex,Ey,num_div)
%SENSORSURFACEGRID Grid on sensor surface
%   Detailed explanation goes here

hlen = 14*1e-3; % sensor size is 28 mm * 28mm
num_ch = size(Rs,2);

grid = linspace(-hlen,hlen,num_div);
num_sensor_grid = num_div.*num_div;
% grid on each sensor coordinate
G = zeros(2,num_sensor_grid);
G(1,:) = reshape(repmat(grid,[num_div,1]),[1,num_sensor_grid]);
G(2,:) = reshape(repmat(grid',[1,num_div]),[1,num_sensor_grid]);

% Grid on sensor surface
Rg = reshape(bsxfun(@plus,[Ex(:),Ey(:)]*G,Rs(:)),[3,num_ch,num_sensor_grid]);
end

