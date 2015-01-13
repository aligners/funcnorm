function spher = computeSphericalFromCartesian(cart, coordMaps)
% FUNCTION computeSphericalFromCartesian
% Returns spher as (phi, theta) pair

% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

N = size(cart, 2);
spher = zeros(2, N, class(cart));
% spher is [phi; theta]

mags = zeros(1, N, class(cart));
try
    mags = sqrt(sum(cart.^2));
catch
    for j = 1:N
        mags(j) = sqrt(sum(cart(:, j).^2));
    end
end

if length(find(imag(cart))) > 0
    % Then something went wrong
    error('computeSphericalFromCartesian:  cart array has imaginary values');
end

xs = cart(1, :);
ys = cart(2, :);
zs = cart(3, :);

if length(coordMaps) == 1
    if coordMaps == 1
        % Then measure phi from the x-axis and theta from the z-axis
        spher(1, :) = acos(xs./mags);
        spher(2, :) = atan2(ys, zs);
    elseif coordMaps == 2
        % Then measure phi from the y-axis and theta from the x-axis
        spher(1, :) = acos(ys./mags);
        spher(2, :) = atan2(zs, xs);
    elseif coordMaps == 3
        % Then measure phi from the z-axis and theta from the x-axis
        spher(1, :) = acos(zs./mags);
        spher(2, :) = atan2(ys, xs);
    else
        error('computeSphericalFromCartesian:  Improper coordinate map');
    end
else
    numFound = 0;
    
    % Find coordinate maps equal to 1 (phi measured from x-axis; theta
    % measured from z-axis)
    inds = find(coordMaps == 1);
    if length(inds) > 0
        numFound = numFound + length(inds);
        spher(1, inds) = acos(xs(inds)./mags(inds));
        spher(2, inds) = atan2(ys(inds), zs(inds));
    end
    
    % Find coordinate map equal to 2 (phi measured from y-axis; theta
    % measured from x-axis
    inds = find(coordMaps == 2);
    if length(inds) > 0
        numFound = numFound + length(inds);
        spher(1, inds) = acos(ys(inds)./mags(inds));
        spher(2, inds) = atan2(zs(inds), xs(inds));
    end
    
    % Find coordinate map equal to 3 (phi measured from z-axis; theta
    % measured from x-axis
    inds = find(coordMaps == 3);
    if length(inds) > 0
        numFound = numFound + length(inds);
        spher(1, inds) = acos(zs(inds)./mags(inds));
        spher(2, inds) = atan2(ys(inds), xs(inds));
    end    
    
    if numFound < N
        error('computeSphericalFromCartesian:  Improper coordinate map');
    end
end