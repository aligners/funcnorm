function cartCoords = normalizeCartesianCoordinates(cartCoords, rho)
% FUNCTION cartCoords = normalizeCartesianCoordinates(cartCoords, rho)
%
% *** INPUT ARGUMENTS ***
%	cartCoords:  a 3 x N matrix of cartesian coordinates
%	rho:  the radius of the sphere to normalize to
%
% *** OUTPUT ARGUMENTS ***
%	cartCoords:  a 3 x N matrix of cartesian coordinates, normalized to a sphere of radius rho
%
%
% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

if nargin == 1
    rho = 1.0;
end

N = size(cartCoords, 2);
try
	cartCoords = cartCoords./repmat(sqrt(sum(cartCoords.^2)), 3, 1);
catch
	for j = 1:N
    	cartCoords(:, j) = rho*cartCoords(:, j)/sqrt(sum(cartCoords(:, j).^2));
	end
end