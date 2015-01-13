function coords = generateMultiHemisphereCoordsDataset(coords, nhem)
% FUNCTION coords = generateMultiHemisphereCoordsDataset(coords, nhem)
%
% *** INPUT ARGUMENTS ***
% coords:  a coordinates structure (see parseSurfaceFile for details)
% nhem:  the number of hemispheres (usually 2)
%
% *** OUTPUT ARGUMENTS ***
% coords: a new coordinates structure for both hemispheres
%
%
% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory


numNodesPerHem = size(coords.cartCoords, 2);
numTrisPerHem = size(coords.triangles, 2);

% First, replicate the cartesian coordinates
coords.cartCoords = repmat(coords.cartCoords, 1, nhem);

% We need to renumber the node numbers of each additional hemisphere
% The original hemisphere is numbered from 0 to (numNodesPerHem - 1)
% The second hemisphere is numbered from (numNodesPerHem) to (2*numNodesPerHem - 1)
for hemNum = 2:nhem
    nbrs2 = coords.neighbors(:, 1:numNodesPerHem);
    % Don't change the -99 nbrs values
    locs = find(nbrs2 ~= -99);
    nbrs2(locs) = nbrs2(locs) + (hemNum-1)*numNodesPerHem;

    coords.neighbors = [coords.neighbors, nbrs2];
    
    tris2 = coords.triangles(:, 1:numTrisPerHem);
    tris2 = tris2 + (hemNum-1)*numNodesPerHem;
    coords.triangles = [coords.triangles, tris2];
end
