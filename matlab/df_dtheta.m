function dVals = df_dtheta(spherCoord1, spherCoord2, resolution, gds)
% FUNCTION dVals = df_dtheta(spherCoord1, spherCoord2, resolution [, gds])
% specify geodesic distances gds if you don't want it re-computed
%
%
% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

if nargin == 3
    % Then we need to compute gds
    gds = computeGeodesicDistances(spherCoord1, spherCoord2);
end

if size(spherCoord1, 2) == 1
	subMode = 'single';
else
	subMode = 'multiple';
	if size(spherCoord1, 2) ~= size(spherCoord2, 2)
		error('df_dtheta:  sphercoord1 must be same size as sphercoord2 in multi submode');
	end
end

h = 0.0201*resolution;

phis1 = spherCoord1(1, :); thetas1 = spherCoord1(2, :);
phis2 = spherCoord2(1, :); thetas2 = spherCoord2(2, :);

% Now any gds that exceed 2*arcsin(h/2) must be set to zero
nonZeroLocs = find(gds < 2*asin(h/2));
dVals = zeros(1, length(gds), class(gds));
if strcmp(subMode, 'single')
	dVals(nonZeroLocs) = -20*(1 - 2*sin(gds(nonZeroLocs)/2)/h).^3.*sin(phis1).*sin(phis2(nonZeroLocs)).*sin(thetas1 - thetas2(nonZeroLocs))/h^2;
else
	dVals(nonZeroLocs) = -20*(1 - 2*sin(gds(nonZeroLocs)/2)/h).^3.*sin(phis1(nonZeroLocs)).*sin(phis2(nonZeroLocs)).*sin(thetas1(nonZeroLocs) - thetas2(nonZeroLocs))/h^2;
end