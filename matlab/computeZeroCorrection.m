function warp_zero = computeZeroCorrection(cartCoords, warps)
% FUNCTION warp_zero = computeZeroCorrection(cartCoords, warps)
% cartCoords should be a 3 x N matrix of cartesian coordinates for each
% node
% warps should be a cell array of warps for each subject
% i.e., warps{1} is a 3 x N matrix

% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

N = size(warps{1}, 2);

if size(cartCoords, 2) == N
    % don't do anything
else
    error('cartCoords does not match size of warp fields');
end

nsubj = length(warps);

avgWarp = zeros(3, N);
for j = 1:nsubj
    warpCartCoords{j} = cartCoords + warps{j};
    avgWarp = avgWarp + warpCartCoords{j}/nsubj;
end

minProj = 1;
for j = 1:N
    % Now project average warp back to the sphere
    avgWarp(:, j) = avgWarp(:, j)/sqrt(sum(avgWarp(:, j).^2));

    % Now we need to compute the rotation matrix that takes avgWarp(:, j)
    % into cartCoords(:, j)
    if avgWarp(:, j)'*cartCoords(:, j) == 1
        % They are too close to rotate
        continue;
    end
    
    % Compute cross product
    c = cross(avgWarp(:, j), cartCoords(:, j));
    c = c/sqrt(sum(c.^2));
    % Now rotate c onto z-axis (rotate about z-axis followed by rotation
    % about y-axis)
    theta = atan2(c(2), c(1));
    phi = acos(c(3));
    R1 = [ cos(theta), sin(theta), 0;
         -sin(theta), cos(theta), 0;
                   0,          0, 1];
    R1 = [ cos(phi), 0, -sin(phi);
                  0, 1,         0;
           sin(phi), 0,  cos(phi)]*R1;

    theta = acos(avgWarp(:, j)'*cartCoords(:, j));
    R2 = [ cos(theta), -sin(theta), 0;
           sin(theta), cos(theta), 0;
                    0,          0, 1];

    R = R1^-1*R2*R1;
    proj = cartCoords(:, j)'*R*avgWarp(:, j);

    % Now use the rotation matrix to rotate each subject
    % This is a really inefficient way to do it...oh well
    for k = 1:nsubj
        w = warpCartCoords{k};
        w(:, j) = R*w(:, j);
        warpCartCoords{k} = w;
    end
end

% Now get warp_zero
for j = 1:nsubj
    warp_zero{j} = warpCartCoords{j} - cartCoords;
end
