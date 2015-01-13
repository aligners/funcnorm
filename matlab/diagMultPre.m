function W = diagMultPre(d, V)
% FUNCTION W = diagMultPre(d, V)
%   Performs W = diag(d)*V

% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

if size(d, 2) ~= 1
    d = d';
end

W = repmat(d, 1, size(V, 2)).*V;
    
