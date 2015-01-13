function W = diagMultPost(V, d)
% FUNCTION W = diagMultPre(V, d)
%   Performs W = V*diag(d)

% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

if size(d, 1) ~= 1
    d = d';
end

W = V.*repmat(d, size(V, 1), 1);
    
