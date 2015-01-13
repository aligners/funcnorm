function spherCoords = computeSphericalCoordsArray(cartCoords, nodeNums, coordMaps)
% This function HASN'T been tested!!!

spherCoords = computeSphericalFromCartesian(cartCoords(:, nodeNums), coordMaps);
%    spherCoordsCellArr{1} = computeSphericalFromCartesian(cartCoords, 1);
%    spherCoordsCellArr{2} = computeSphericalFromCartesian(cartCoords, 2);
%    spherCoordsCellArr{3} = computeSphericalFromCartesian(cartCoords, 3);


%% warp cartcoords
%wcc = reshape(repmat(warpCartCoords, maxNumNbrs, 1), 3, N*maxNumNbrs);
%wnc = warpCartCoords(:, fullUpdNbrs(:)+1);
