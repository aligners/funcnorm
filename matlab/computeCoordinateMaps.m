function coordMaps = computeCoordinateMaps(cartCoords, warpCartCoords)
% The objective of this function is to find an appropriate coordinate
% mapping of the sphere for each voxel such that the voxel's original
% coordinate location and current warped location are not located close to
% the poles (phi = 0, phi = pi)

% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

N = size(cartCoords, 2);

coordMaps = zeros(1, N, 'uint8');

for j = 1:N
    % For each voxel, 
    c1 = max(abs(cartCoords(1, j)), abs(warpCartCoords(1, j)));
    c2 = max(abs(cartCoords(2, j)), abs(warpCartCoords(2, j)));
    c3 = max(abs(cartCoords(3, j)), abs(warpCartCoords(3, j)));

    [garbage, coordMaps(j)] = min([c1 c2 c3]);
end
