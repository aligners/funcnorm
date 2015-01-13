function gds = computeGeodesicDistances(coord1, coord2, rho)
% FUNCTION gds = computeGeodesicDistances(coord1, coords2, rho)
%   coord1 are the coordinates for one voxel
%   If coord1 is a 3 x 1 vector, then coordinates are Cartesian (x, y, z)
%   If coord1 is a 2 x 1 vector, then coordinates are spherical (phi, theta)
%   coord2 are the coordinates of N voxels 
%   coord2 must  be in the same coordinate representation as coord1, meaning
%   size(coord1, 1) == size(coord2, 1)
%   So coord2 is either 2 x N or 3 x N, depending on coord1
%
%   This function returns a 1 x N vector representing the spherical geodesic
%   distance of coord1 to each of the voxels in coord2
%
% *** OR ***
%
% FUNCTION gds = computeGeodesicDistances(coord1, coord2, rho)
%   coord1 are the coordinates for N voxels
%   If coord1 is 3 x N vector, then coordinates are Cartesian (x, y, z)
%   If coord1 is 2 x N vector, then coordinates are spherical (phi, theta)
%   coord2 are the coordinates of N voxels as well
%   coord2 must be in the same coordinate representation as coord1, meaning
%   size(coord1, 1) == size(coord2, 1)
%   
%   This function returns a 1 x N vector representing the spherical
%   geodesic distances between the columns of coord1 with the cooresponding
%   column in coord2
%
%
% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

if nargin == 2
    rho = 1;
end

if size(coord1, 1) == 2
    mode = 'spherical';
elseif size(coord1, 1) == 3
    mode = 'cartesian';
else
    error('computeGeodesicDistances:  coord1 must be 2 x N or 3 x N');
end

if size(coord1, 2) == 1
    subMode = 'single';
else
    subMode = 'multi';
    if size(coord1, 2) ~= size(coord2, 2)
        error('computeGeodesicDistances:  coord1 and coord2 must have the same number of columns more multi subMode');
    end
end

if size(coord1, 1) ~= size(coord2, 1)
    error('computeGeodesicDistances:  coord1 and coord2 must be in the same coordinate representation');
end

N = size(coord2, 2);
% Allocate memory for return argument
gds = zeros(1, N, class(coord1));

if strcmp(mode, 'spherical')
    % We are in the spherical mode
    if strcmp(subMode, 'single')
        % Get the (phi, theta) spherical coordinates from coord1
        phi1 = coord1(1); theta1 = coord1(2);
        % Get the (phi, theta) spherical coordinates for all of the voxels in
        % coord2
        phis2 = coord2(1, :); thetas2 = coord2(2, :);
    
        % A straightforward approach to computing geodesic distance would be
        % through the arccos of the dot product
        % However, arccos is numerically unstable for small geodesic distances
        % So we use this formula instead, which utilizes arcsin instead of
        % arccos
        gds = 2*rho*asin(sqrt(sin(phi1)*sin(phis2).*(sin((theta1 - thetas2)/2)).^2 + (sin((phi1 - phis2)/2)).^2));
        % If we were not on the unit sphere (but on a sphere with radius rho), then gds would be computed as:
        % gds = 2*rho*asin(sqrt(sin(phi1)*sin(phis2).*(sin((theta1 - thetas2)/2)).^2 + (sin((phi1 - phis2)/2)).^2));
    else
        % Get the (phi, theta) spherical coordinates from coord1
        phis1 = coord1(1, :); thetas1 = coord1(2, :);
        % Get the (phi, theta) spherical coordinates for all of the voxels in
        % coord2
        phis2 = coord2(1, :); thetas2 = coord2(2, :);
    
        % A straightforward approach to computing geodesic distance would be
        % through the arccos of the dot product
        % However, arccos is numerically unstable for small geodesic distances
        % So we use this formula instead, which utilizes arcsin instead of
        % arccos
        gds = 2*rho*asin(sqrt(sin(phis1).*sin(phis2).*(sin((thetas1 - thetas2)/2)).^2 + (sin((phis1 - phis2)/2)).^2));
        % If we were not on the unit sphere (but on a sphere with radius rho), then gds would be computed as:
        % gds = 2*rho*asin(sqrt(sin(phi1)*sin(phis2).*(sin((theta1 - thetas2)/2)).^2 + (sin((phi1 - phis2)/2)).^2));
    end
else
    % We are in the cartesian mode
    if strcmp(subMode, 'single')
        % A straightforward approach to computing geodesic distances would
        % be through the arccos of the dot product
        % However, arccos is numerically unstable for small geodesic
        % distances
        % So we use this formula instead, which utilizes arcsin instead of
        % arccos
        % This uses the arcsin of the magnitude of the cross product
        gds = rho*asin(sqrt(sum(cross(repmat(coord1, 1, N), coord2).^2)));
    else
        gds = rho*asin(sqrt(sum(cross(coord1, coord2).^2)));
    end
end
