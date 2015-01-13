function spherWarp = convertCartesianWarpToSphericalWarp(cartCoords, cartWarp, coordMaps)

% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

N = size(cartCoords, 2);
spherWarp = zeros(2, N, class(cartCoords));

for j = 1:N
    currSpherCoords = computeSphericalFromCartesian(cartCoords(:, j), coordMaps(:, j));
    currWarpSpherCoords = computeSphericalFromCartesian(cartCoords(:, j) + cartWarp(:, j), coordMaps(:, j));

    spherWarp(:, j) = currWarpSpherCoords - currSpherCoords;
    % Should check to make sure that there are no out of bounds changes
end
