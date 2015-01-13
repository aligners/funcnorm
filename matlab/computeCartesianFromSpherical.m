function cart = computeCartesianFromSpherical(spher, coordMaps)

% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

N = size(spher, 2);
cart = zeros(3, N, class(spher));

if length(find(imag(spher))) > 0
    error('computeCartesianFromSpherical:  spher array contains imaginary components');
end

phis = spher(1, :);
thetas = spher(2, :);

if length(coordMaps) == 1
    % Then use the same coordinate map for all
    coordMaps = coordMaps*ones(1, N);

    if coordMaps == 1
        % Then phi is measured from x-axis and theta from z-axis
        cart(1, :) = cos(phis);
        cart(2, :) = sin(phis).*sin(thetas);
        cart(3, :) = sin(phis).*cos(thetas);
    elseif coordMaps == 2
        % Then phi is measured from the y-axis and theta from the x-axis
        cart(1, :) = sin(phis).*cos(thetas);
        cart(2, :) = cos(phis);
        cart(3, :) = sin(phis).*sin(thetas);
    elseif coordMaps == 3
        % Then phis is measured from the z-axis and theta from the x-axis
        cart(1, :) = sin(phis).*cos(thetas);
        cart(2, :) = sin(phis).*sin(thetas);
        cart(3, :) = cos(phis);
    else
        error('computeCartesianFromSpherical:  Improper coordinate map');
    end
else
    numFound = 0;
    
    % Do coordMaps == 1
    inds = find(coordMaps == 1);
    if length(inds) > 0
        numFound = numFound + length(inds);
        cart(1, inds) = cos(phis(inds));
        cart(2, inds) = sin(phis(inds)).*sin(thetas(inds));
        cart(3, inds) = sin(phis(inds)).*cos(thetas(inds));
    end
    
    % Do coordMaps == 2
    inds = find(coordMaps == 2);
    if length(inds) > 0
        numFound = numFound + length(inds);
        cart(1, inds) = sin(phis(inds)).*cos(thetas(inds));
        cart(2, inds) = cos(phis(inds));
        cart(3, inds) = sin(phis(inds)).*sin(thetas(inds));
    end
    
    % Do coordMaps == 3
    inds = find(coordMaps == 3);
    if length(inds) > 0
        numFound = numFound + length(inds);
        cart(1, inds) = sin(phis(inds)).*cos(thetas(inds));
        cart(2, inds) = sin(phis(inds)).*sin(thetas(inds));
        cart(3, inds) = cos(phis(inds));
    end
    
    if numFound < N
        error('computeCartesianFromSpherical:  Improper coordinate map');
    end
end
