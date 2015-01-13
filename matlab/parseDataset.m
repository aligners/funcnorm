function [U, sigma, V, M, N] = parseDataset(dataset)

% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory
  
U = [];
if strcmp(class(dataset), 'struct')
    % Then we already have all the requisite matrices
    if isfield(dataset, 'sigma')
        sigma = dataset.sigma;
    else
        error('Dataset structure does not contain a sigma field');
    end

    if isfield(dataset, 'V')
        V = dataset.V;
    else
        error('Dataset structure does not contain a V field');
    end

    % Now compute M and N
    [m, n] = size(sigma);
    % We want sigma to be a column vector
    if m == 1
        sigma = sigma';
        M = n;
    elseif n == 1
        M = m;
    else
        error('Dataset sigma must be a column vector');
    end

    [n, m] = size(V);
    % Make sure it's of the correct size
    if m ~= M
        error('Dataset V field must contain the same number of columns as the length of the sigma field');
    end
    N = n;

    if isfield(dataset, 'U') 
        U = dataset.U;
        if size(U, 1) ~= M || size(U, 2) ~= M
            error('Dataset U field must be a square matrix with the same number of rows as the length of the sigma field');
        end
    end 
else
    % Otherwise, we have to compute the SVD ourselves
    [M, N] = size(dataset);
    [U, S, V] = svd(dataset, 'econ');
    sigma = diag(S);
end
    