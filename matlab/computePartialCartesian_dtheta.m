function dp_dtheta = computePartialCartesian_dtheta(cart, coordMaps)
% FUNCTION dp_dtheta = computePartialCartesian_dtheta(cart, coordMaps)
% dp is a 3 x N matrix, where N = size(cart, 2)

% Copyright (c) 2008 Bryan Conroy, Jim Haxby, Peter Ramadge, Mert Rory Sabuncu, and Ben Singer

N = size(cart, 2);

if length(coordMaps) == 1
    coordMaps = repmat(coordMaps, 1, N);
end

dp_dtheta = zeros(3, N, class(cart));

for j = 1:N
    x = cart(1, j);
    y = cart(2, j);
    z = cart(3, j);
    mag = sqrt(x^2 + y^2 + z^2);

    if coordMaps(j) == 1
        phi = acos(x/mag);
        theta = atan2(y, z);

        dp_dtheta(1, j) = 0;
        dp_dtheta(2, j) = sin(phi)*cos(theta);
        dp_dtheta(3, j) = -sin(phi)*sin(theta);
    elseif coordMaps(j) == 2
        phi = acos(y/mag);
        theta = atan2(z, x);

        dp_dtheta(1, j) = -sin(phi)*sin(theta);
        dp_dtheta(2, j) = 0;
        dp_dtheta(3, j) = sin(phi)*cos(theta);
    elseif coordMaps(j) == 3
        phi = acos(z/mag);
        theta = atan2(y, x);

        dp_dtheta(1, j) = -sin(phi)*sin(theta);
        dp_dtheta(2, j) = sin(phi)*cos(theta);
        dp_dtheta(3, j) = 0;
    end
end
