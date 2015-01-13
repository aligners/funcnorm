function [triArea, dareal_dphi, dareal_dtheta] = computeArealTerms(triangles, cartCoords, coordMaps, origTriAreas, triNormals, rho)
% FUNCTION computeArealTerms

% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

if nargin == 5
    rho = 1;
end

N = size(cartCoords, 2);
T = size(triangles, 2);

if  nargout == 1
    computeDerivatives = 0;
    triArea = cast(0, class(cartCoords));
    dareal_dphi = [];
    dareal_dtheta = [];
else
    computeDerivatives = 1;
    triArea = cast(0, class(cartCoords));
    dareal_dphi = zeros(N, 1, class(cartCoords));
    dareal_dtheta = zeros(N, 1, class(cartCoords));
    dp_dphi = computePartialCartesian_dphi(cartCoords, coordMaps);
    dp_dtheta = computePartialCartesian_dtheta(cartCoords, coordMaps);
end

triIndices = triangles+1;
p0s = rho*cartCoords(:, triIndices(1, :));
p1s = rho*cartCoords(:, triIndices(2, :));
p2s = rho*cartCoords(:, triIndices(3, :));
  
as = p1s - p0s;
bs = p2s - p0s;
  
abcrosses = cross(as, bs);
  
newAreas = sum(abcrosses.*triNormals)/2;
% Find locations where newAreas is less than zero (indicates folding of the cortex)
locs = find(newAreas <= 0);
aTerms = transpose(newAreas - origTriAreas);
triArea = sum(aTerms(locs).^2);
 
if computeDerivatives == 1
	bcns = cross(bs, triNormals);
 	ncas = cross(triNormals, as);
 	
 	for l = 1:length(locs)
 		loc = locs(l);
 		
 		cI = triangles(1, loc) + 1;
		dp0_dphi = dp_dphi(:, cI);
		dp0_dtheta = dp_dtheta(:, cI);
		dareal_dphi(cI) = dareal_dphi(cI) - aTerms(loc)*(bcns(:, loc) + ncas(:, loc))'*dp0_dphi;
		dareal_dtheta(cI) = dareal_dtheta(cI) - aTerms(loc)*(bcns(:, loc) + ncas(:, loc))'*dp0_dtheta;
		
		cI = triangles(2, loc) + 1;
		dp1_dphi = dp_dphi(:, cI);
		dp1_dtheta = dp_dtheta(:, cI);
		dareal_dphi(cI) = dareal_dphi(cI) + aTerms(loc)*bcns(:, loc)'*dp1_dphi;
		dareal_dtheta(cI) = dareal_dtheta(cI) + aTerms(loc)*bcns(:, loc)'*dp1_dtheta;
		
		cI = triangles(3, loc) + 1;
		dp2_dphi = dp_dphi(:, cI);
		dp2_dtheta = dp_dtheta(:, cI);
		dareal_dphi(cI) = dareal_dphi(cI) + aTerms(loc)*ncas(:, loc)'*dp2_dphi;
		dareal_dtheta(cI) = dareal_dtheta(cI) + aTerms(loc)*ncas(:, loc)'*dp2_dtheta;
	end
end

%% OLD, non-vectorized method
%for j = 1:T
%  currNodes = triangles(:, j);
%
%  % First compute the triarea
%  p0 = rho*cartCoords(:, currNodes(1)+1);
%  p1 = rho*cartCoords(:, currNodes(2)+1);
%  p2 = rho*cartCoords(:, currNodes(3)+1);
%
%  a = p1 - p0;
%  b = p2 - p0;
%  n = triNormals(:, j);
%
%  newArea = cross(a, b)'*n/2;
%
%  % First, if the triangle area is positive, then just keep going
%  if newArea > 0
%      continue;
%  end
%
%  aTerm = newArea - origTriAreas(j);
%  triArea = triArea + aTerm^2;
%
%  if computeDerivatives == 0
%     continue;
%  end
%
%  bcn = cross(b, n);
%  nca = cross(n, a);
%
%  % Do node 1 first
%  cI = currNodes(1)+1;
%  dp0_dphi = dp_dphi(:, cI);
%  dp0_dtheta = dp_dtheta(:, cI);
%  dareal_dphi(cI) = dareal_dphi(cI) - aTerm*(bcn + nca)'*dp0_dphi;
%  dareal_dtheta(cI) = dareal_dtheta(cI) - aTerm*(bcn + nca)'*dp0_dtheta;
%
%  % Now do node 2
%  cI = currNodes(2)+1;
%  dp1_dphi = dp_dphi(:, cI);
%  dp1_dtheta = dp_dtheta(:, cI);
%  dareal_dphi(cI) = dareal_dphi(cI) + aTerm*bcn'*dp1_dphi;
%  dareal_dtheta(cI) = dareal_dtheta(cI) + aTerm*bcn'*dp1_dtheta;
%
%  % Now do node 3
%  cI = currNodes(3)+1;
%  dp2_dphi = dp_dphi(:, cI);
%  dp2_dtheta = dp_dtheta(:, cI);
%  dareal_dphi(cI) = dareal_dphi(cI) + aTerm*nca'*dp2_dphi;
%  dareal_dtheta(cI) = dareal_dtheta(cI) + aTerm*nca'*dp2_dtheta;
%end
%% END OLD, non-vectorized method
