function dp_dphi = computePartialCartesian_dphi(cart, coordMaps)
% FUNCTION dp_dphi = computePartialCartesian_dphi(cart, coordMaps)
% dp is a 3 x N matrix, where N = size(cart, 2)
    
% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

N = size(cart, 2);

if length(coordMaps) == 1
    coordMaps = repmat(coordMaps, 1, N);
end

dp_dphi = zeros(3, N, class(cart));

for j = 1:N
    x = cart(1, j);
    y = cart(2, j);
    z = cart(3, j);
    mag = sqrt(x^2 + y^2 + z^2);

    if coordMaps(j) == 1
        phi = acos(x/mag);
        theta = atan2(y, z);

        dp_dphi(1, j) = -sin(phi);
        dp_dphi(2, j) = cos(phi)*sin(theta);
        dp_dphi(3, j) = cos(phi)*cos(theta);
    elseif coordMaps(j) == 2
        phi = acos(y/mag);
        theta = atan2(z, x);

        dp_dphi(1, j) = cos(phi)*cos(theta);
        dp_dphi(2, j) = -sin(phi);
        dp_dphi(3, j) = cos(phi)*sin(theta);
    elseif coordMaps(j) == 3
        phi = acos(z/mag);
        theta = atan2(y, x);

        dp_dphi(1, j) = cos(phi)*cos(theta);
        dp_dphi(2, j) = cos(phi)*sin(theta);
        dp_dphi(3, j) = -sin(phi);
    end
end
