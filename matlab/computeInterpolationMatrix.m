function [A, Alocs, locsLength, gds, Q, dAD_dphi, dAD_dtheta] = computeInterpolationMatrix
	
% updNbrs is (maxNumNbrs x numNodes)
maxNumNbrs = size(updNbrs, 1);
N = size(updNbrs, 2);
opsDataType = 'single';

if nargout == 1
    computeDerivative = 0;
    dAD_dphi = [];
    dAD_dtheta = [];
elseif nargout == 3
    computeDerivative = 1;
    dAD_dphi = -99*ones(maxNumNbrs, N, opsDataType);
    dAD_dtheta = -99*ones(maxNumNbrs, N, opsDataType);
else
    error('computeInterpolationMatrix:  Improper call');
end

A = sparse(N, N);
Alocs = -99*ones(maxNbrs, N, 'int32');
locsLength = zeros(N, 1, 'int32');

% Find the locations where no neighbor exists (i.e., updNbrs == -99)
neg99Locs = find(updNbrs == -99);
fullUpdNbrs = updNbrs;
fullUpdNbrs(neg99Locs) = floor(neg99Locs/maxNumNbrs) - 1;

% warp cartcoords
wcc = reshape(repmat(warpCartCoords, maxNumNbrs, 1), 3, N*maxNumNbrs);
% warp neighbor cartcoords
wnc = warpCartCoords(:, fullUpdNbrs(:)+1);

gds = computeGeodesicDistances(wcc, wnc);
AVals = interp_f(wcc, wnc, res, gds, opsDataType);
AVals(neg99Locs) = 0;

jCopy = reshape(repmat(1:N, maxNumNbrs, 1), 1, N*maxNumNbrs);
A(fullUpdNbrs(:) + 1, jCopy) = AVals;

Q = V1ST*A;
V1SBig = V1S(Alocs, 1:);


% Now normalize Q
try
    D = norm(Q);
catch
	D = zeros(1, N);
	for j = 1:N
		D(j) = norm(Q(:, j));
	end
end

% Now we want to normalize the columns of Q -- perform post-multiplication
Q = diagMultPost(Q, D);

if computeDerivative == 1
	cM = reshape(repmat(coordMaps, maxNumNbrs, 1), 1, N*maxNumNbrs);
	spherCoords1 = computeSphericalFromCartesian(wcc, cM);
	spherCoords2 = computeSphericalFromCartesian(wnc, cM);

    % Now we must compute gradient information as well
    % This computes derivative of interpolation coefficients with
    % respect to phi and theta
    dA_dphi = transpose(df_dphi(currWarpSpherCoords, nbrSpherCoords, res, gds));
    dA_dtheta = transpose(df_dtheta(currWarpSpherCoords, nbrSpherCoords, res, gds));
            
            Q_V1S = Q(:, j)'*V1ST(:, NI);
            % This computes the derivative of the normalizing weight
            % (weight which normalizes data matrix to unit norm) with
            % respect to both phi and theta
            dD_dphi = D^2*Q_V1S*dA_dphi;
            dD_dtheta = D^2*Q_V1S*dA_dtheta;
            % This computes derivative of A*D (normalized interpolation
            % matrix) with respect to phi and theta
            dAD_dphi(1:currLength, j) = D*dA_dphi - dD_dphi*A;
            dAD_dtheta(1:currLength, j) = D*dA_dtheta - dD_dtheta*A;                
        end

    
    
        % Compute geodesic distance to all neighbors
        gds = computeGeodesicDistances(currWarpSpherCoords, nbrSpherCoords);
        % Now compute interpolation function
        [A, nonZeroLocs] = interp_f(currWarpSpherCoords, nbrSpherCoords, res, gds, opsDataType);
        
        % Now we must trim everything based on nonZeroLocs
        % We also want A to be a column vector
        A = A(nonZeroLocs)';
        currNbrs = currNbrs(nonZeroLocs);
        nbrSpherCoords = nbrSpherCoords(:, nonZeroLocs);
        gds = gds(nonZeroLocs);
        % Set Alocs (which specifies rows where each column of A is
        % nonzero)
        currLength = length(currNbrs);
        locsLength(j) = currLength;
        NI = currNbrs+1;
        Alocs(1:currLength, j) = NI;

        Q(:, j) = V1ST(:, NI)*A;
        % Normalize the column to unit norm
        qnorm = norm(Q(:, j));
        if qnorm <= 1e-10
            D = 0;
        else
            D = qnorm^-1;
        end
        Q(:, j) = D*Q(:, j);

        if nargout >= 2
            % Now we must compute gradient information as well
            % This computes derivative of interpolation coefficients with
            % respect to phi and theta
            dA_dphi = transpose(df_dphi(currWarpSpherCoords, nbrSpherCoords, res, gds));
            dA_dtheta = transpose(df_dtheta(currWarpSpherCoords, nbrSpherCoords, res, gds));
            
            Q_V1S = Q(:, j)'*V1ST(:, NI);
            % This computes the derivative of the normalizing weight
            % (weight which normalizes data matrix to unit norm) with
            % respect to both phi and theta
            dD_dphi = D^2*Q_V1S*dA_dphi;
            dD_dtheta = D^2*Q_V1S*dA_dtheta;
            % This computes derivative of A*D (normalized interpolation
            % matrix) with respect to phi and theta
            dAD_dphi(1:currLength, j) = D*dA_dphi - dD_dphi*A;
            dAD_dtheta(1:currLength, j) = D*dA_dtheta - dD_dtheta*A;                
        end
    end
	displayLogItem(sprintf('\t\t%s', 'Completed computing interpolation matrix A!'), logFile);
