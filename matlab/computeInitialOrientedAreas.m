function [triAreas, orientedNormals] = computeInitialOrientedAreas(triangles, cartCoords, rho)
% FUNCTION [triAreas, orientedNormals] = computeInitialOrientedAreas(triangles, cartCoords, rho)
% cartCoords should be unit norm, but specify rho to use a different
% spherical radius

% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

if nargin == 2
    rho = 1.0;
end

numTriangles = size(triangles, 2);

triIndices = triangles+1;

p0s = cartCoords(:, triIndices(1, :));
p1s = cartCoords(:, triIndices(2, :));
p2s = cartCoords(:, triIndices(3, :));

as = rho*(p1s - p0s);
bs = rho*(p2s - p0s);

ns = cross(as, bs);


nNorms = sqrt(sum(ns.^2));
triAreas = nNorms/2;

try
	% Try the vectorized method first
	orientedNormals = ns./repmat(nNorms, 3, 1);
catch
	% Sometimes the vectorized method fails for out of memory problems
	orientedNormals = zeros(3, numTriangles, class(ns));
	for j = 1:size(ns, 2)
		orientedNormals(:, j) = ns(:, j)/sqrt(sum(ns(:, j).^2));
	end
end


%% OLD CODE -- non-vectorized method
% triAreas = zeros(1, numTriangles, class(cartCoords));
% orientedNormals = zeros(3, numTriangles, class(cartCoords));
% 
% % Now compute the oriented areas for each triangle
% for j = 1:numTriangles
%    % First, compute the oriented normal
%    % p0 is the main vertex of the triangle
%    p0 = cartCoords(:, triangles(1, j)+1);
%    p1 = cartCoords(:, triangles(2, j)+1);
%    p2 = cartCoords(:, triangles(3, j)+1);
% 
%    a = rho*(p1 - p0);
%    b = rho*(p2 - p0);
% 
%    n = cross(a, b);
% 
%    nNorm = norm(n, 2);
%    orientedNormals(:, j) = n/nNorm;
%    triAreas(j) = nNorm/2;
% end
%% END OLD CODE -- non-vectorized method