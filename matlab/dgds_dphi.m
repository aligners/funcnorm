function dVals = dgds_dphi(spherCoords1, spherCoords2, resolution, gds)
% FUNCTION dVals = dgds_dphi(spherCoord1, spherCoords2, resolution, gds)
%
% *** OR ***
%
% FUNCTION dVals = dgds_dphi(spherCoords1, spherCoords2, resolution, gds)
%
%
% Compute derivative of geodesic distance between spherCoords1 and
% spherCoords2 with respect to phi spherical coordinates of spherCoords1
%
% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

if size(spherCoords1, 2) == 1
    spherCoords1 = repmat(spherCoords1, 1, size(spherCoords2, 2));
end

if nargin == 3
    % Then we need to compute gds
    gds = computeGeodesicDistances(spherCoord1, spherCoords2);
end

h = 0.0201*resolution;

phis1 = spherCoords1(1, :); thetas1 = spherCoords1(2, :);
phis2 = spherCoords2(1, :); thetas2 = spherCoords2(2, :);

% Now any gds that are too small should be kept at zero
nonZeroLocs = find(gds > 1e-5);
dVals = zeros(1, length(gds), class(gds));
dVals(nonZeroLocs) = (sin(phis1(nonZeroLocs)).*cos(phis2(nonZeroLocs)) - cos(phis1(nonZeroLocs)).*sin(phis2(nonZeroLocs)).*cos(thetas1(nonZeroLocs) - thetas2(nonZeroLocs)))./sin(gds(nonZeroLocs));
