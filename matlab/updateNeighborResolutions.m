function [retNbrs, retResNbrSizes, retTotalNbrs] = updateNeighborResolutions(cartCoords, warpCartCoords, nbrs, resNbrSizes, totalNbrs, updNbrs, updResNbrSizes, updTotalNbrs)
% FUNCTION [retNbrs, retResNbrSizes, retTotalNbrs] = updateNeighborResolutions(cartCoords, warpCartCoords, nbrs, resNbrSizes, totalNbrs, updNbrs, updResNbrSizes, updTotalNbrs)
% This function is compatible with nodes that lie on multiple hemispheres

% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

N = size(cartCoords, 2);

% Allocate matrices for the updated neighbor information
retNbrs = -99*ones(size(nbrs), 'int32');
retResNbrSizes = zeros(size(resNbrSizes), 'int32');
retTotalNbrs = zeros(size(totalNbrs), 'int32');

for j = 1:N
    % First, initialize the ret nbrs to the current nbrs
    try
        retNbrs(:, j) = updNbrs(:, j);
    catch
        [j, size(retNbrs), size(updNbrs)]
        rethrow(lasterror);
    end
    retResNbrSizes(:, j) = updResNbrSizes(:, j);
    retTotalNbrs(j) = updTotalNbrs(j);

    closestNbrPre = retNbrs(1, j);
    while 1
        % Get the current (at this iteration) neighbors for this node
        currNbrs = retNbrs(1:retTotalNbrs(j), j);
        nbrsCartCoords = cartCoords(:, currNbrs+1);

        % Want to find the closest neighbor grid point -- find the one
        % with maximal projection
        projections = warpCartCoords(:, j)'*nbrsCartCoords;
        [maxVal, I] = max(projections);

        % Now use the original neighbors data structure to assign the
        % updated neighbors
        retNbrs(:, j) = nbrs(:, currNbrs(I)+1);
        retResNbrSizes(:, j) = resNbrSizes(:, currNbrs(I)+1);
        retTotalNbrs(j) = totalNbrs(currNbrs(I)+1);

        % Now check if we need to do another iteration
        % If closest neighbor number hasn't changed from the past
        % iteration, then we are done
        if closestNbrPre == currNbrs(I)
            break;
        else
            closestNbrPre = currNbrs(I);
        end
    end
end