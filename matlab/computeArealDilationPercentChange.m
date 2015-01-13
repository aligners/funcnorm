function pctChange = computeArealDilationPercentChange(cartCoords, warpCartCoords, triangles, rho)
% FUNCTION pctChange = computeArealDilationPercentChange(cartCoords, warpCartCoords, triangles)

% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

N = size(cartCoords, 2);
T = size(triangles, 2);

triIndices = triangles+1;

% Compute old areas
p0s = cartCoords(:, triIndices(1, :));
p1s = cartCoords(:, triIndices(2, :));
p2s = cartCoords(:, triIndices(3, :));
as = p1s - p0s;
bs = p2s - p0s;
abcrosses = cross(as, bs);
abNorms = sqrt(sum(abcrosses.^2));
oldAreas = abNorms/2;

% Compute new areas
p0s = warpCartCoords(:, triIndices(1, :));
p1s = warpCartCoords(:, triIndices(2, :));
p2s = warpCartCoords(:, triIndices(3, :));
as = p1s - p0s;
bs = p2s - p0s;
abcrosses = cross(as, bs);
abNorms = sqrt(sum(abcrosses.^2));
newAreas = abNorms/2;

oldNodeAreas = zeros(1, N);
newNodeAreas = zeros(1, N);
for j = 1:T
	oldNodeAreas(triIndices(:, j)) = oldNodeAreas(triIndices(:, j)) + oldAreas(j);
	newNodeAreas(triIndices(:, j)) = newNodeAreas(triIndices(:, j)) + newAreas(j);
end

sum(oldNodeAreas)
sum(newNodeAreas)

pctChange = 100*(newNodeAreas - oldNodeAreas)./oldNodeAreas;
