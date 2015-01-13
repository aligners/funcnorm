function TW = computeInterpOnSphere(T, cartCoords, neighbors, warp, nn)
% nn specifies nearest-neighbor interpolation
% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

N = size(T, 2);
if ~exist('nn'); nn=0; end;

warpCartCoords = cartCoords + warp;
TW = zeros(size(T));

numNbrs = zeros(1, N);
neighbors = [0:(N-1); neighbors];
for j = 1:N
    numNbrs(j) = length(find(neighbors(:, j) ~= -99));
end

avgNb = 0;
for j = 1:N
    currCartCoords = warpCartCoords(:, j);
    
    % The first step is to find the closest grid point
    % This will be done in an iterative fashion
    % Initialize 
    cNbrs = neighbors(1:numNbrs(j), j);
    prevNbr = j-1;
    
    while 1
        projs = currCartCoords'*cartCoords(:, cNbrs+1);
        [maxVals, I] = max(projs);
        closestNeighbor = cNbrs(I);
        
        cNbrs = neighbors(1:numNbrs(closestNeighbor+1), closestNeighbor+1);
        
        if prevNbr == closestNeighbor
            break;
        else
            prevNbr = closestNeighbor;
        end
    end
    
    % Now we can find the other two closest grid points among the current
    % neighbors
    [maxVals, I] = sort(currCartCoords'*cartCoords(:, cNbrs+1), 2, 'descend');
    
    if nn == 1
    	triNbrs = cNbrs(I(1));
    else
	    triNbrs = cNbrs(I(1:3));
	end
    
    triCartCoords = cartCoords(:, triNbrs+1);
    gds = computeGeodesicDistances(currCartCoords, triCartCoords);
    [interpVals, nonZeroLocs] = interp_f(currCartCoords, triCartCoords, 1, gds);
    if length(nonZeroLocs) == 0
        error(['problem here:  ', num2str(j)]);
    end
    avgNb = avgNb + length(nonZeroLocs);
    
    TW(:, j) = T(:, triNbrs+1)*interpVals';
end
end