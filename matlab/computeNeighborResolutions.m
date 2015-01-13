function [nbrs, resNbrSizes, totalNbrs] = computeNeighborResolutions(nearestNbrs, maxRes)
% FUNCTION computeNeighborResolutions(nearestNbrs, maxRes)
% This function is compatible with nodes that lie on multiple hemispheres

% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

[nbrsRows, N] = size(nearestNbrs);

% Initially allocate a large matrix for neighbors (we'll trim it down
% later)
nbrs = -99*ones(200, N, 'int32');
% resNbrSizes specifies the number of neighbors added in each resolution
resNbrSizes = zeros(maxRes, N, 'int32');
% totalNeighbors specifies the total number of neighbors (over all
% resolutions) for each voxel
totalNbrs = zeros(1, N, 'int32');

% This data structure (resPtrs) specifies the start and end indices for the
% neighbors added in the previous resolution
resPtrs = zeros(2, N, 'int32');
% First, work on resolution #1
for j = 1:N
    % Account for itself being a neighbor
    s = length(find(nearestNbrs(:, j) ~= -99));
    if nearestNbrs(1, j) ~= (j-1)
        s = s + 1;
        nbrs(1:s, j) = [(j-1); nearestNbrs(1:(s-1), j)];
    else
        nbrs(1:s, j) = nearestNbrs(1:s, j);
    end

    resNbrSizes(1, j) = s;
    resPtrs(:, j) = [1; s];
    totalNbrs(j) = totalNbrs(j) + s;
end

% Now, iterate over larger resolutions
for res = 2:maxRes
    % Within each resolution, iterate over voxels
    for j = 1:N
        prevStart = resPtrs(1, j); prevEnd = resPtrs(2, j);
        % Account for the fact that the first element is always itself
        if prevStart == 1
            prevStart = 2;
        end

        % Get the neighbors added from the previous resolution
        prevNbrs = nbrs(prevStart:prevEnd, j);

        % Get the neighbors of these neighbors
        currNbrs = unique(reshape(nearestNbrs(:, prevNbrs+1), length(prevNbrs)*nbrsRows, 1));
        % Remove any dummy values
        currNbrs(currNbrs == -99) = [];
        % Now we only want to keep the new neighbors
        % Keep the values in currNbrs that are not in nbrs(1:prevEnd, j)
        % nbrs(1:prevEnd, j) contains all non-dummy neighbors of the
        % current voxel
        currNbrs = setdiff(currNbrs, nbrs(1:prevEnd, j));

        % Update the size for this resolution
        newSize = length(currNbrs);
        resNbrSizes(res, j) = newSize;
        % Also update the overallSize
        totalNbrs(j) = totalNbrs(j) + newSize;

        % We also need to update the resPtrs for the next time around
        newStart = prevEnd + 1; newEnd = newStart + newSize - 1;
        resPtrs(:, j) = [newStart; newEnd];

        % Lastly, update the resolution neighbors matrix
        nbrs(newStart:newEnd, j) = currNbrs;
    end
end

% Now trim the neighbors matrix
maxNbrs = max(totalNbrs');
% Trim the neighbors matrix down to the appropriate size
nbrs = nbrs(1:maxNbrs, :);