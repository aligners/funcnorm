function warp = funcnorm_register(warpDataset, templateDataset, coords, x, regularization, maxResolution, logFile)
% FUNCTION warp = funcnorm_register(warpDataset, templateDataset, coords, x, regularization, maxResolution [, logFile])
%
% Let M = time-series dimension
% Let N = # cortical nodes
%
% *** INPUT ARGUMENTS ***
% warpDataset:  the dataset to be warped
% templateDataset:  the template dataset
%       These datasets should be M x N matrices with node time-series as
%       columns
% coords:  A structure specifying the (initial pre-warp) spatial coordinates of the cortical nodes
%       The structure has the following fields:
%           i)   .cartCoords:  Cartesian coordinates (3 x N matrix), in which
%                      each column specifies an (x, y, z) coordinate
%           ii)  .neighbors:  neighbors of each voxel (6 x N matrix),
%                      in which each column specifies the node numbers of the 6
%                      closest neighbors to a voxel
%           iii) .triangles:  matrix (3 x N) of node neighbors
% x: initial warp field
%       This should be specified as a 3 x N matrix, with the jth column
%       specifying the (x, y, z) warp for the jth cortical node
% regularization:  structure which specifies how to regularize the warp field
%       Examples of the structure:
%       a) Based on metric and areal distortion:
%           i)   .mode = 'metric_and_areal_distortion'
%           ii)  .lambda_metric = 1 (weighting of metric distortion term in objective function)
%           iii) .lambda_areal = 1 (weighting of areal distortion term in objective function)
%       b) No regularization:
%           i) .mode = 'none'
% maxResolution:  the starting resolution of the algorithm (using multi-resolution approach)
% logFile (optional):  filename for log items
%
% *** OUTPUT ARGUMENT ***
% warp:  3 x N matrix specifying the optimized warp field for the
%   templateDataset

% Copyright (c) 2008 Bryan Conroy, Jim Haxby, Peter Ramadge, Mert Rory Sabuncu, and Ben Singer

    % Set the radius of the sphere (we'll work on unit sphere)
    rho = 1;
    
    if nargin == 6
        logFile = '';
    end    
    
%% Parse the datasets
    displayLogItem('Parsing the datasets...', logFile);
    % Variable suffix of 1 corresponds to the warp dataset (W1, C1, U1, s1, S1, V1)
    % Variable suffix of 2 corresponds to the template dataset (W2, C2, U2, s2, S2, V2)
    % Parse the datasets
    % **** Warp dataset ****
    % Let W1 = timeseries matrix (M x N) of the warp dataset
    % Then W1 = U1*diag(sigma1)*V1'
    %         = U1*S1*V1'
    % **** Template dataset ****
    % Let W2 = timeseries matrix (M x N) of the template dataset
    % Then W2 = U2*diag(sigma2)*V2'
    %         = U2*S2*V2'
    [U1, s1, V1, M1, N1] = parseDataset(warpDataset);
    clear warpDataset;
    % Make sure all new variables are cast in the original datatype
    opsDataType = class(U1);
    [U2, s2, V2, M2, N2] = parseDataset(templateDataset);
    clear templateDataset;
    if N1 ~= N2
        error('Number of cortical nodes must match between the two datasets');
    end
    N = N1;
    clear N1 N2;
    
    if M1 ~= M2
        error('Number of time-series points must match between the two datasets');
    end
    M = max(length(s1), length(s2));
    clear M1 M2;
    
    V1S = diagMultPost(V1, s1); % V1S = V1*diag(s1)
    clear s1 V1;
    
    V2S = diagMultPost(V2, s2); % V2S = V2*diag(s2)

    V2Orig = V2;
    s2Orig = s2;
    U2Orig = U2;
    displayLogItem('Completed parsing the datasets', logFile);
%% End parsing the datasets


%% Parse the coords
    displayLogItem('Parsing the coordinates and neighbors structures...', logFile);
    % Now we have all the information in the structure
    if isfield(coords, 'neighbors')
        nbrs = coords.neighbors;
    else
        error('coords structure does not contain a neighbors field');
    end
        
    hasTriangles = 0;
    if isfield(coords, 'triangles')
        triangles = coords.triangles;
        hasTriangles = 1;
    end
	% Don't error here because it may not be necessary to specify the
	% triangles field
        
	if isfield(coords, 'cartCoords')
        cartCoords = normalizeCartesianCoordinates(coords.cartCoords, rho);
        if size(cartCoords, 1) ~= 3 || size(cartCoords, 2) ~= N
            error(sprintf('Improper cartCoords field specified in coords (must be of size 3 x N, N=%d)', N));
        end                
    else
        error('The coords structure does not contain a cartCoords field');
    end         
    clear coords;
    
    % Now compute the spherical coordinates based in all coordinate maps
    spherCoordsCellArr{1} = computeSphericalFromCartesian(cartCoords, 1);
    spherCoordsCellArr{2} = computeSphericalFromCartesian(cartCoords, 2);
    spherCoordsCellArr{3} = computeSphericalFromCartesian(cartCoords, 3);
        
    % Now compute the neighbors at all resolutions
    origNbrs = nbrs;
    origNumNbrs = zeros(1, N, 'int32');
    for j = 1:N
        origNumNbrs(j) = length(find(origNbrs(:, j) ~= -99));
    end
    [nbrs, resNbrSizes, totalNbrs] = computeNeighborResolutions(nbrs, maxResolution);
    maxNbrs = max(totalNbrs);
    displayLogItem('Completed parsing the coordinates and neighbors structures', logFile);
%% End parsing the coords


%% REGULARIZATION
    if strcmp(regularization.mode, 'metric_and_areal_distortion')
        displayLogItem('Using metric and areal distortion terms for regularization...', logFile);
        if isfield(regularization, 'lambda_metric') == 0
            error('regularization structure must contain lambda_metric field');
        end
        if isfield(regularization, 'lambda_areal') == 0
            error('regularization structure must contain lambda_areal field');
        end

        % We need to compute the initial metric distances
        if regularization.lambda_metric > 0
            regularization.lambda_metric = regularization.lambda_metric/(4*N);
            
            regularization.metricDistances = -99*ones(6, N, opsDataType);
            displayLogItem('Computing initial metric distances...', logFile);
            for j = 1:N
                currCartCoords = cartCoords(:, j);
                numNbrs = origNumNbrs(j);            
                currNbrs = origNbrs(1:numNbrs, j);
                nbrCartCoords = cartCoords(:, currNbrs+1);
            
                regularization.metricDistances(1:numNbrs, j) = transpose(computeGeodesicDistances(currCartCoords, nbrCartCoords));
            end
            displayLogItem('Completed computing initial metric distances!', logFile);
        end
        
        if regularization.lambda_areal > 0
            % Make sure that the triangles field has been specified
            if hasTriangles == 0
                error('In order to compute areal distortion terms, the triangles must be specified in the coords structure');
            end

            % Now we must compute the initial oriented area terms
            displayLogItem('Computing initial oriented area terms...', logFile);
            regularization.triangles = triangles;
            clear triangles hasTriangles;

            numTriangles = size(regularization.triangles, 2);
            regularization.lambda_areal = regularization.lambda_areal/(2*numTriangles);
            [regularization.orientedAreas, regularization.orientedNormals] = computeInitialOrientedAreas(regularization.triangles, cartCoords);
            displayLogItem('Completed computing initial oriented area terms', logFile);
        end
    elseif strcmp(regularization.mode, 'none')
        displayLogItem('No regularization used', logFile);
    else
        error(['Regularization mode: ''', regularization.mode, ''' is not recognized.']);
    end
%% END REGULARIZATION


    optMethod = 'steepest_descent';    
    displayLogItem(['Using ', optMethod, ' as optimization technique'], logFile);

    % Initialize the updated neighbor information
    updNbrs = nbrs;
    updResNbrSizes = resNbrSizes;
    updTotalNbrs = totalNbrs;
    
    % Before we begin, we must convert the 3 x N warp matrix x to a set of
    % (dphi, dtheta) warp parameters
    warpCartCoords = cartCoords + x;
    coordMaps = computeCoordinateMaps(cartCoords, warpCartCoords);
    x = convertCartesianWarpToSphericalWarp(cartCoords, x, coordMaps);
    x = x(:);

    
%% Main Loop for Optimization
% Calls fminunc at each iteration
    for res = maxResolution:-1:1
        displayLogItem(['Beginning pass at resolution #', num2str(res)], logFile);
        Q = V1S';
        V1ST = V1S';
        
        if res >= 4
            tolFun = 1e-1;
        elseif res == 3
            tolFun = 1e-2;
        elseif res == 2
            tolFun = 1e-3;
        elseif res == 1
            tolFun = 1e-3;
        end
        
        % Use this for steepest descent
        maxIter = 300;
        maxFunEvals = 500;
        options = optimset('MaxIter', maxIter, 'MaxFunEvals', maxFunEvals, 'TolFun', tolFun, 'OutputFcn', @outfn, 'DerivativeCheck', 'off', 'GradObj', 'on', 'Display', 'iter', 'LargeScale', 'on', 'Hessian', 'on');
        displayLogItem(['Parameters used in this optimization:  MaxIter = ', num2str(maxIter), '; MaxFunEvals = ', num2str(maxFunEvals), '; TolFun = ', num2str(tolFun)], logFile);

        displayLogItem(['Blurring the template dataset at resolution #', num2str(res)], logFile);
        % Then we must blur the template dataset
        [V2, s2, U2] = blurDataset(V2Orig, s2Orig, U2Orig, cartCoords, nbrs, totalNbrs, res);
        displayLogItem(['Completed blurring the template dataset at resolution #', num2str(res)], logFile);
        
        W2TU1 = diagMultPost(V2, s2)*(U2'*U1);

        % Reset the warped coordinate locations
        x = reshape(x, 2, N);
        warpCartCoords = computeWarpCoordsFromSphericalWarp(cartCoords, x, coordMaps);
        % We need to determine the correct coordinate mapping for each
        % node
        coordMaps = computeCoordinateMaps(cartCoords, warpCartCoords);
        x = x(:);
        
        % fminunc only works with data types of class double
        [x, f, exitflag, output] = fminunc(@computeObjective, cast(x, 'double'), options);
        displayLogItem(['Results of optimization at resolution #', num2str(res), ':'], logFile);
            displayLogItem(sprintf('\t%s', ['Number of iterations:  ', num2str(output.iterations)]), logFile);            
            displayLogItem(sprintf('\t%s', ['Number of function evaluations:  ', num2str(output.funcCount)]), logFile);
            displayLogItem(sprintf('\t%s', ['Algorithm used:  ', output.algorithm]), logFile);
            displayLogItem(sprintf('\t%s', ['Exit message:  ', output.message]), logFile);
        
        % Now we have to trim neighbors from the outermost resolution
        % First, decrease the neighbors count for each voxel by the number
        % of voxel neighbors in the outermost resolution
        totalNbrs = totalNbrs - resNbrSizes(end, :);
        % Next, remove the neighbor count for the outermost resolution
        resNbrSizes = resNbrSizes(1:end-1, :);
        
        updTotalNbrs = updTotalNbrs - updResNbrSizes(end, :);
        updResNbrSizes = updResNbrSizes(1:end-1, :);
        
        % We could also re-trim the neighbors matrix
        maxNbrs = max(max(totalNbrs), max(updTotalNbrs));
        nbrs = nbrs(1:maxNbrs, :);
        updNbrs = updNbrs(1:maxNbrs, :);
    end
%% END Main Loop for Optimization
    
    x = reshape(x, 2, N);
    warpCartCoords = computeWarpCoordsFromSphericalWarp(cartCoords, x, coordMaps);
    warp = warpCartCoords - cartCoords;
    
    return;

                                                                                    
                                                                                    
%% ************************************************************************
%   FUNCTION stop = outfn(x, optimVals, state)
    function stop = outfn(x, optimVals, state)        
        stop = false;
        if strcmp(state, 'iter') == 0 && strcmp(state, 'init') == 0
            return;
        end

        displayLogItem(['Iteration number:  ', num2str(optimVals.iteration)], logFile);
        
        corrAvg = 0;
        for j = 1:N
            corrAvg = corrAvg + W2TU1(j, :)*Q(:, j)/N;
        end

        displayLogItem(sprintf('\t%s', ['Average inter-subject correlation:  ', num2str(corrAvg)]), logFile);        
    end;
%% END FUNCTION outfn                                                                                    
        

%% ************************************************************************
%   FUNCTION [f,g,H] = computeObjective(x)
%   x specifies the current (dphi, dtheta) warp for each voxel
function [f,g,H] = computeObjective(x)
    % This function is a nested function of connectivity_align
    % It uses the following global variables:
    %   nbrs, resNbrSizes, totalNbrs, maxNbrs
    %   cartCoords
    
    x = cast(reshape(x, 2, N), opsDataType);
	pause(0.01);

    % First, compute the updated warped cartesian coordinates
    warpCartCoords = computeWarpCoordsFromSphericalWarp(cartCoords, x, coordMaps);
    % We need to determine the correct coordinate mapping for each node
    coordMaps = computeCoordinateMaps(cartCoords, warpCartCoords);
    % Update the neighbors data structure 
    [updNbrs, updResNbrSizes, updTotalNbrs] = updateNeighborResolutions(cartCoords, warpCartCoords, nbrs, resNbrSizes, totalNbrs, updNbrs, updResNbrSizes, updTotalNbrs);
    % Now compute the warped locations in spherical coordinates
    warpSpherCoords = computeSphericalFromCartesian(warpCartCoords, coordMaps);
    
	% Re-initialize all data structures
	Q = zeros(M, N, opsDataType);
	if nargout >= 2
        % Then we must compute gradient -- allocate space for it
        dAD_dphi = zeros(maxNbrs, N, opsDataType);
        dAD_dtheta = zeros(maxNbrs, N, opsDataType);
        
        g = zeros(2*N, 1, opsDataType);
    end
	if nargout >= 3
        % Then we must compute Hessian -- allocate space for it
        H = sparse(2*N, 2*N);
        H = spdiags(ones(2*N, 1), 0, H);
    end
        
	Alocs = -99*ones(maxNbrs, N, 'int32');
	locsLength = zeros(N, 1, 'int32');
        
	% Compute A matrix
	% If nargout >= 2, then also compute dA_dx and dA_dy
	% For each j, compute the jth row of A
	for j = 1:N
        currCoordMap = coordMaps(:, j);
        spherCoords = spherCoordsCellArr{currCoordMap};
        
        currWarpSpherCoords = warpSpherCoords(:, j);

        % Retrieve the current neighbors
        currNbrs = updNbrs(1:updTotalNbrs(j), j);
        % Get the spherical coordinates of all of the neighbors
        nbrSpherCoords = spherCoords(:, currNbrs+1);
        
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
        
    f = 0;        
    mult_intersubj = 1.0;%/N;
    for j = 1:N
        f = f + mult_intersubj*(1-W2TU1(j, :)*Q(:, j));
	end
    
    if strcmp(regularization.mode, 'metric_and_areal_distortion')
        % We need to compute the areal and metric distortion
        if regularization.lambda_metric > 0
            % First, compute metric distortion terms
            count = 1;
            for j = 1:N
                currWarpCartCoords = warpCartCoords(:, j);
                % Get its original 6 nearest neighbors
                currNumNbrs = origNumNbrs(j);
                currNbrs = origNbrs(1:currNumNbrs, j);
                nbrWarpCartCoords = warpCartCoords(:, currNbrs+1);
            
                currInitMD = regularization.metricDistances(1:currNumNbrs, j);

                % We need to get the spherical coordinates with respect to
                % the current coordinate mapping
                currCoordMap = coordMaps(:, j);
                spherCoords = spherCoordsCellArr{currCoordMap};
                currWarpSpherCoords = spherCoords(:, j) + x(:, j);
                nbrWarpSpherCoords = computeSphericalFromCartesian(nbrWarpCartCoords, currCoordMap);
                gds = computeGeodesicDistances(currWarpSpherCoords, nbrWarpSpherCoords);
                currMD = gds';
                mdDiffs = currMD - currInitMD;
            
                f = f + regularization.lambda_metric*sum(mdDiffs.^2);
            
                if nargout >= 2                
                    % Compute the gradient for this term as well here
                    dg_dphi = dgds_dphi(currWarpSpherCoords, nbrWarpSpherCoords, res, gds);
                    % dg_dphi is a row vector, while mdDiffs is a column vector
                    g(count) = g(count) + 4*regularization.lambda_metric*(dg_dphi*mdDiffs);
                    dg_dtheta = dgds_dtheta(currWarpSpherCoords, nbrWarpSpherCoords, res, gds);
                    g(count+1) = g(count + 1) + 4*regularization.lambda_metric*(dg_dtheta*mdDiffs);
                    count = count + 2;
                end
            end
        end
        
        if regularization.lambda_areal > 0
            % Now compute the areal term
            if nargout >= 2
                [triArea, dareal_dphi, dareal_dtheta] = computeArealTerms(regularization.triangles, warpCartCoords, coordMaps, regularization.orientedAreas, regularization.orientedNormals);
                g(1:2:end) = g(1:2:end) + 2*regularization.lambda_areal*dareal_dphi;
                g(2:2:end) = g(2:2:end) + 2*regularization.lambda_areal*dareal_dtheta;
            else
                triArea = computeArealTerms(regularization.triangles, warpCartCoords, coordMaps, regularization.orientedAreas, regularization.orientedNormals);                
            end
            f = f + regularization.lambda_areal*triArea;
        end
    end

    if nargout >= 2
        % Then we need to compute the gradient
        count = 1;
        for p = 1:N
            locLength = locsLength(p);
            locs = Alocs(1:locLength, p);

            dS_dAD = W2TU1(p, :)*V1S(locs, :)';
            g(count) = g(count) - mult_intersubj*dS_dAD*dAD_dphi(1:locLength, p);
            g(count + 1) = g(count + 1) - mult_intersubj*dS_dAD*dAD_dtheta(1:locLength, p);

            count = count + 2;
        end
    end
    if isa(f,'double') == 0
        f = double(f);
    end
    if isa(g,'double') == 0
        g = double(g);
    end
    g = real(g);
    f = real(f);
end
%% ** END FUNCTION computeObjective
end
%% ** END funcnorm_register
