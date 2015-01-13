function [interpVals, nonZeroLocs] = interp_f(coord1, coord2, resolution, gds, opsDataType, zM)
% FUNCTION [interpVals, nonZeroLocs] = interp_f(coord1, coords2, resolution [[, gds] [, opsDataType] [, zM]])
% coord1, coord2 can be specified as cartesian or spherical coordinates
% returns a row vector interpVals as well as a vector nonZeroLocs which
% specifies locations where the interpolation values are nonzero
% specify geodesic distances gds if you don't want it re-computed
% opsDataType is optional
% zM is an optional zero-mask (indices to ignore in coord2)
%
% *** OR ***
%
% FUNCTION [interpVals, nonZeroLocs] = interp_f(coord1, coord2, resolution, gds, opsDataType)
%	coord1 are the coordinates for N voxels
%	If coord1 is 3 x N vector, then coordinates are Cartesian (x,y,z)
%	If coord1 is 2 x N vector, then coordinates are spherical (phi, theta)
%	coord2 are the coordinates of N voxels as well
%	coord2 must be in the same coordinate representation as coord1, meaning size(coord1, 1) == size(coord2, 1)
%	specify gds if you don't want it re-computed
%
%	This function returns a 1 x N vector (interpVals) representing the interpolation coefficients of coord2 centered at coord1
%	This function also returns a vector nonZeroLocs, which specifies the locations in interpVals that are nonzero
%
%
% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

if size(coord1, 1) == 2
    mode = 'spherical';
elseif size(coord1, 1) == 3
    mode = 'cartesian';
end

if size(coord1, 2) == 1
    subMode = 'single';
else
    subMode = 'multi';
    if size(coord1, 2) ~= size(coord2, 2)
        error('interp_f:  coord1 and coord2 must have the same number of columns more multi subMode');
    end
end


if nargin == 3
    % Then we need to compute gds
    gds = computeGeodesicDistances(coord1, coord2);
end

if nargin < 5
    opsDataType = 'double';
end

if nargin < 6
	zM = [];
end
	
% average gds distance to closest neighbors is 0.0201
h = 0.0201*resolution;

% Now any gds that exceed 2*arcsin(h/2) must be set to zero
if length(zM) > 0	
	nonZeroLocs = setdiff(find(gds < 2*asin(h/2)), zM);
else
	nonZeroLocs = find(gds < 2*asin(h/2));
end

interpVals = zeros(1, length(gds), opsDataType);
s = 2*sin(gds(nonZeroLocs)/2)/h;
interpVals(nonZeroLocs) = (1 - s).^4.*(4*s + 1);