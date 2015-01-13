function warpCartCoords = computeWarpCoordsFromSphericalWarp(cartCoords, spherWarp, coordMaps)

% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

N = size(cartCoords, 2);

spherCoords = computeSphericalFromCartesian(cartCoords, coordMaps);
warpSpherCoords = spherCoords + spherWarp;
warpCartCoords = computeCartesianFromSpherical(warpSpherCoords, coordMaps);

%warpCartCoords = zeros(3, N, class(cartCoords));
%for j = 1:N
%    currSpherCoords = computeSphericalFromCartesian(cartCoords(:, j), coordMaps(:, j));
%    currWarpSpherCoords = currSpherCoords + spherWarp(:, j);
%    warpCartCoords(:, j) = computeCartesianFromSpherical(currWarpSpherCoords, coordMaps(:, j));
%end
