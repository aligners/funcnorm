function [mdDiff, dmdDiffs_dphi, dmdDiffs_dtheta] = computeMetricTerms(origNbrs, origNumNbrs, cartCoords, coordMaps, origMetricDistances, res, rho)
% FUNCTION computeMetricTerms

% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

N = size(cartCoords, 2);
maxNumNbrs = size(origNbrs, 1);

if nargout == 1
    computeDerivatives = 0;
    dmdDiffs_dphi = [];
    dmdDiffs_dtheta = [];
elseif nargout == 3
    computeDerivatives = 1;
    dmdDiffs_dphi = zeros(N, 1, class(cartCoords));
    dmdDiffs_dtheta = zeros(N, 1, class(cartCoords));
else
    error('computeMetricTerms:  Number of output arguments can be 1 or 3');
end

if nargin == 6
    rho = 1;
elseif nargin == 7
    % do nothing
else
    error('computeMetricTerms:  improper call');
end

locs = find(origNbrs == -99);
fullOrigNbrs = origNbrs;
fullOrigNbrs(locs) = floor((locs-1)/maxNumNbrs); % OLD:  NOT RIGHT --> floor(locs/maxNumNbrs) - 1;
cc = reshape(repmat(cartCoords, maxNumNbrs, 1), 3, N*maxNumNbrs);
nc = cartCoords(:, fullOrigNbrs(:)+1);

% The first column of currMD corresponds to the distances of node 1 to its
% neighbors
gds = computeGeodesicDistances(cc, nc, rho);
currMD = reshape(gds, maxNumNbrs, N);

mdDiffsMat = currMD - origMetricDistances;
mdDiff = sum(sum(mdDiffsMat.^2));

if computeDerivatives == 1
    cM = reshape(repmat(coordMaps, maxNumNbrs, 1), 1, N*maxNumNbrs);
    
    csc = computeSphericalFromCartesian(cc, cM);
    nsc = computeSphericalFromCartesian(nc, cM);

    dg_dphi = reshape(dgds_dphi(csc, nsc, res, gds), maxNumNbrs, N);
    dg_dtheta = reshape(dgds_dtheta(csc, nsc, res, gds), maxNumNbrs, N);

    dmdDiffs_dphi = transpose(dot(dg_dphi, mdDiffsMat));
    dmdDiffs_dtheta = transpose(dot(dg_dtheta, mdDiffsMat));
end


%% BEGIN OLD CODE:  non-vectorized version%%
% for j = 1:N
%     currWarpCartCoords = cartCoords(:, j);
%     currNumNbrs = origNumNbrs(j);
%     currNbrs = origNbrs(1:currNumNbrs, j);
%     nbrWarpCartCoords = cartCoords(:, currNbrs+1);
%     
%     currInitMD = origMetricDistances(1:currNumNbrs, j);
%     
%     currCoordMap = coordMaps(:, j);
%     currWarpSpherCoords = computeSphericalFromCartesian(currWarpCartCoords, currCoordMap);
%     nbrWarpSpherCoords = computeSphericalFromCartesian(nbrWarpCartCoords, currCoordMap);
%     gds = computeGeodesicDistances(currWarpSpherCoords, nbrWarpSpherCoords);
%     currMD = gds';
%     mdDiffsMat = currMD - currInitMD;
%     
%     mdDiffs(j) = sum(mdDiffsMat.^2);
%     
%     dg_dphi = dgds_dphi(currWarpSpherCoords, nbrWarpSpherCoords, res, gds);
%     dmdDiffs_dphi(j) = dg_dphi*mdDiffsMat;
%     
%     dg_dtheta = dgds_dtheta(currWarpSpherCoords, nbrWarpSpherCoords, res, gds);
%     dmdDiffs_dtheta(j) = dg_dtheta*mdDiffsMat;
% end
%% END OLD CODE:  non-vectorized version %%