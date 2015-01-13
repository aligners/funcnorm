function [V, s, U] = blurDataset(V, s, varargin)
% FUNCTION [V, s, U] = blurDataset(V, s, U, cartCoords, nbrs, totalNbrs, resolution)
% ** OR **
% FUNCTION [V, s] = blurDataset(V, s, cartCoords, nbrs, totalNbrs, resolution)

% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

if length(varargin) == 5
    if nargout == 3
        returnU = 1;
        U = varargin{1};
    else
        returnU = 0;
    end
    
    cartCoords = varargin{2};
    nbrs = varargin{3};
    totalNbrs = varargin{4};
    resolution = varargin{5};
elseif length(varargin) == 4
    if nargout == 3
        error('blurTemplate:  Cannot return U matrix if it is not supplied in the input arguments');
    end

    returnU = 0;
    cartCoords = varargin{1};
    nbrs = varargin{2};
    totalNbrs = varargin{3};
    resolution = varargin{4};
else
    error('blurTemplate: Improper number of input arguments.  Type ''help blurTemplate'' for function prototype');
end

[N, M] = size(V);

Y = diagMultPre(s, V'); % Y = diag(s)*V'
opsDataType = class(Y);


% BEGIN:  New interpolation code
maxNbrs = cast(max(totalNbrs), 'double');
Q = zeros(M, N, opsDataType);
% Find the locations where no neighbor exists (i.e., updNbrs == -99)
neg99Locs = find(nbrs == -99);
fullUpdNbrs = nbrs;
fullUpdNbrs(neg99Locs) = floor((neg99Locs-1)/maxNbrs);
cc = reshape(repmat(cartCoords, maxNbrs, 1), 3, N*maxNbrs);
nc = cartCoords(:, fullUpdNbrs(:)+1);
% Compute geodesic distances between coordinates
gds = computeGeodesicDistances(cc, nc);
% Avals will be N*maxNbrs x 1
[AVals, nonZeroLocs] = interp_f(cc, nc, resolution, gds, opsDataType, neg99Locs);
NI = fullUpdNbrs(nonZeroLocs); % NI are the nonzero neighbor indices
A = reshape(AVals, maxNbrs, N);        
ACell = mat2cell(A, maxNbrs, ones(1, N));
% Q(:, j) = V1ST(:, NI)*A;
NSplitSize = ceil(N/maxNbrs);
nStart = NSplitSize*(0:(maxNbrs-1));
nStart(nStart > N) = [];
for ns = nStart
	ne = min(ns + NSplitSize, N);
    NC = (maxNbrs*ns:(maxNbrs*ne - 1)) + 1;
	YCell = mat2cell(Y(1:M, fullUpdNbrs(NC)+1), M, maxNbrs*ones(1, (ne-ns)));
    Q(1:M, (ns+1):ne) = cellMultPost(YCell, ACell((ns+1):ne));
    mags = sqrt(sum(Q(:, (ns+1):ne).^2));
    mLocs = find(mags > 1e-8);
    D = zeros(1, (ne-ns));
    D(mLocs) = 1./mags(mLocs);
    Q(1:M, (ns+1):ne) = diagMultPost(Q(:, (ns+1):ne), D);
end
% END:  New interpolation code    


% BEGIN:  OLD INTERPOLATION CODE
%Q = zeros(M, N, opsDataType);
%for j = 1:N
%     currCartCoords = cartCoords(:, j);
%     numNbrs = totalNbrs(j);
%     currNbrs = nbrs(1:numNbrs, j);
%     nbrCartCoords = cartCoords(:, currNbrs+1);
% 
%     % Compute geodesic distances on sphere to all neighboring nodes
%     % Compute interpolation function at neighbor node locations
%     [A, nonZeroLocs] = interp_f(currCartCoords, nbrCartCoords, resolution);
%     A = A(nonZeroLocs)';
%     currNbrs = currNbrs(nonZeroLocs);
%     Q(:, j) = Y(:, currNbrs+1)*A;
%     qnorm = norm(Q(:, j));
%     if qnorm > 1e-10
%         Q(:, j) = Q(:, j)/qnorm;
%     end
% end
% END:  OLD INTERPOLATION CODE




% Now compute SVD of Q matrix
[Unew, S, V] = svd(Q, 'econ');
s = diag(S);

if returnU
    U = U*Unew;
end
