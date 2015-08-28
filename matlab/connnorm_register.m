function warp = connnorm_register(warpDataset, templateDataset, coords, x, regularization, maxResolution, logFile)
% FUNCTION warp = connnorm_register(warpDataset, templateDataset, coords, x, regularization, maxResolution [, logFile])
%
% Let M1, M2 = time-series dimension of warpDataset, templateDataset, respectively
% Let N = # cortical nodes
%
% *** INPUT ARGUMENTS ***
% warpDataset:  the dataset to be warped
% templateDataset:  the template dataset
%       These datasets should be M1 x N and M2 x N matrices with node time-series as
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
%
%
% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

    % Set the radius of the sphere (we'll work on unit sphere)
    rho = 1;

	compMethod = 'fast'; % Either slow or fast
    % The fast method used to  be experimental, but now is verified to be
    % bug-free.  So use the fast method (unless you want the script to take
    % longer!)
        
    if nargin == 6
        logFile = '';
    end



%% Parse the datasets
displayLogItem(['Parsing the datasets...'], logFile);
    % Variable suffix of 1 corresponds to the warp dataset (W1, C1, U1, s1, S1, V1)
    % Variable suffix of 2 corresponds to the template dataset (W2, C2, U2, s2, S2, V2)
    % Parse the datasets
    % **** Warp dataset ****
    % Let W1 = timeseries matrix (M1 x N) of the warp dataset
    % Let C1 = correlation matrix (N x N) of the warp dataset
    % Then W1 = U1*diag(sigma1)*V1'
    %         = U1*S1*V1'
    % Then C1 = V1*diag(sigma1.^2)*V1'
    %         = V1*(S1.^2)*V1'
    % **** Template dataset ****
    % Let W2 = timeseries matrix (M2 x N) of the template dataset
    % Let C2 = correlation matrix (N x N) of the template dataset
    % Then W2 = U2*diag(sigma2)*V2'
    %         = U2*S2*V2'
    % Then C2 = V2*diag(sigma2.^2)*V2'
    %         = V2*(S2.^2)*V2'
    hasU = 0;
    [U1, s1, V1, M1, N1] = parseDataset(warpDataset);
    opsDataType = class(V1);
    if length(U1) > 0
        hasU = 1;
    end
    clear warpDataset;
    [U2, s2, V2, M2, N2] = parseDataset(templateDataset);
    clear templateDataset;
    if N1 ~= N2
        error('Number of cortical nodes must match between the two datasets');
    end
    N = N1;
    clear N1 N2;
    
    msMatch = 1;
    if M1 ~= M2
        msMatch = 0;
        % We can also clear U1 and U2 variables (since they are irrelevant)
        hasU = 0;
        U1 = [];
        U2 = [];
    end

    V1S = diagMultPost(V1, s1); % V1S = V1*diag(s1)
    V1S0 = [V1S; zeros(1, M1)]; % V1S with a zero-row appended
    clear s1 V1;
    
    V2S = diagMultPost(V2, s2); % V2S = V2*diag(s2)

    V2Orig = V2;
    s2Orig = s2;
    U2Orig = U2;
    
    displayLogItem(sprintf('\t%s', ['Number of cortical nodes:  ', num2str(N)]), logFile);            
    displayLogItem(sprintf('\t%s', ['Number of timepoints in warp dataset:  ', num2str(M1)]), logFile);
    displayLogItem(sprintf('\t%s', ['Number of timepoints in template dataset:  ', num2str(M2)]), logFile);
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

    opsCoordsType = class(cartCoords); % Format for coords structures
    
    % Now compute the spherical coordinates based in all coordinate maps
    % Don't need this anymore -- use sC instead (see a few lines below)
%    spherCoordsCellArr{1} = computeSphericalFromCartesian(cartCoords, 1);
%    spherCoordsCellArr{2} = computeSphericalFromCartesian(cartCoords, 2);
%    spherCoordsCellArr{3} = computeSphericalFromCartesian(cartCoords, 3);
    
	sC = zeros(6, N, opsCoordsType);
	sC(1:2, 1:N) = computeSphericalFromCartesian(cartCoords, 1);
	sC(3:4, 1:N) = computeSphericalFromCartesian(cartCoords, 2);
	sC(5:6, 1:N) = computeSphericalFromCartesian(cartCoords, 3);
    
    % Now compute the neighbors at all resolutions
    origNbrs = nbrs;
    origNumNbrs = zeros(1, N, 'int32');
    for j = 1:N
        origNumNbrs(j) = length(find(origNbrs(:, j) ~= -99));
    end
    [nbrs, resNbrSizes, totalNbrs] = computeNeighborResolutions(nbrs, maxResolution);
    maxNbrs = cast(max(totalNbrs), 'double');
displayLogItem('Completed parsing the coordinates and neighbors structures', logFile);
%% End parsing the coords


%% REGULARIZATION
displayLogItem('Parsing the regularization structure...', logFile);
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
            
            regularization.metricDistances = zeros(6, N, opsCoordsType);
            displayLogItem('Computing initial metric distances...', logFile);

            % BEGIN:  NEW CODE TO COMPUTE INITIAL METRIC DISTANCES
            neg99Locs = find(origNbrs == -99);
            fullOrigNbrs = origNbrs;
            fullOrigNbrs(neg99Locs) = floor((neg99Locs-1)/maxNbrs);
            cc = reshape(repmat(cartCoords, size(fullOrigNbrs, 1), 1), 3, N*size(fullOrigNbrs, 1));
            nc = cartCoords(:, fullOrigNbrs(:)+1);

            regularization.metricDistances = reshape(computeGeodesicDistances(cc, nc), size(fullOrigNbrs, 1), N);
            % END:  NEW CODE TO COMPUTE INITIAL METRIC DISTANCES
            
            % BEGIN:  OLD CODE TO COMPUTE INITIAL METRIC DISTANCES
%             for j = 1:N
%                 currCartCoords = cartCoords(:, j);
%                 numNbrs = origNumNbrs(j);            
%                 currNbrs = origNbrs(1:numNbrs, j);
%                 nbrCartCoords = cartCoords(:, currNbrs+1);
%             
%                 regularization.metricDistances(1:numNbrs, j) = transpose(computeGeodesicDistances(currCartCoords, nbrCartCoords));
%             end
            % END:  OLD CODE TO COMPUTE INITIAL METRIC DISTANCES
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
displayLogItem('Completed parsing the regularization structure!', logFile);
%% END REGULARIZATION


    optMethod = 'steepest_descent';    
    displayLogItem(['Using ', optMethod, ' as optimization technique'], logFile);
	lambda_intrasubj = 1.0/N;
	lambda_intersubj = 0.0;
	
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

    % This is for debugging purposes
    % Will check the derivative accuracy
    derivCheck = 'off';
    if strcmp(derivCheck, 'on')
    	displayLogItem(['*** WARNING:  We are in derivative check mode!!!! ***'], logFile);
	    numDerivCheckNodes = 5;
		% Generate (numDerivCheckNodes) random nodes whose derivatives we will check
    	% derivNodeNums ranges between 0 and N-1
		derivNodeNums = round((N-1)*rand(numDerivCheckNodes, 1)); 
		derivWarpInds = zeros(1, 2*numDerivCheckNodes);
		derivWarpInds(1:2:end) = 2*derivNodeNums + 1;
		derivWarpInds(2:2:end) = derivWarpInds(1:2:end) + 1;
		
	    derivWarp = x(derivWarpInds);
    end
    
%% Main Loop for Optimization
% Calls fminunc at each iteration
    for res = maxResolution:-1:1
        displayLogItem(['Beginning pass at resolution #', num2str(res)], logFile);
        Q = V1S';
%        V1ST = V1S';
        
        if res >= 4
            tolFun = 1e-1;
        elseif res == 3
            tolFun = 1e-1;
        elseif res == 2
            tolFun = 1e-2;
        elseif res == 1
            tolFun = 1e-3;
        end
        
        % Use this for steepest descent
        maxIter = 300;
        maxFunEvals = 500;
        % Default is using LargeScale=on when objective function supplies gradient
        % Set 'LargeScale','GradObj', and 'Hessian' off if gradient computation fails. 
        % Ref: fminnunc doc in Matlab
        options = optimset('MaxIter', maxIter, 'MaxFunEvals', maxFunEvals, 'TolFun', tolFun, 'OutputFcn', @outfn, 'DerivativeCheck', derivCheck, 'GradObj', 'on', 'Display', 'iter', 'LargeScale', 'on', 'Hessian', 'on', 'FunValCheck', 'on');
        displayLogItem(['Parameters used in this optimization:  MaxIter = ', num2str(maxIter), '; MaxFunEvals = ', num2str(maxFunEvals), '; TolFun = ', num2str(tolFun)], logFile);

        % Then we must blur the template dataset
        if res > 1
%        	displayLogItem(['WARNING:  NOT blurring the template dataset at resolution #', num2str(res)], logFile);
%        	V2 = V2Orig;
%        	s2 = s2Orig;
	        displayLogItem(['Blurring the template dataset at resolution #', num2str(res)], logFile);
	        if hasU
		        [V2, s2, U2] = blurDataset(V2Orig, s2Orig, U2Orig, cartCoords, nbrs, totalNbrs, res);
	    	    W2TU1 = diagMultPost(V2, s2)*(U2'*U1);
		    else
		    	[V2, s2] = blurDataset(V2Orig, s2Orig, cartCoords, nbrs, totalNbrs, res);
	  		end
	        displayLogItem(['Completed blurring the template dataset at resolution #', num2str(res)], logFile);
	  	else
	  		if hasU
	  			V2 = V2Orig;
	  			s2 = s2Orig;
	  			U2 = U2Orig;
	  			W2TU1 = diagMultPost(V2, s2)*(U2'*U1);
	  		else
	  			V2 = V2Orig;
	  			s2 = s2Orig;
	  		end
	  	end
        
        s2_squared = s2.^2;
        V2S2T = diagMultPre(s2_squared, V2');
        
        if strcmp(derivCheck, 'on')       
        	x = zeros(2*N, 1);
        	x(derivWarpInds) = derivWarp;
        	x = reshape(x, 2, N);
        	
        	warpCartCoords = computeWarpCoordsFromSphericalWarp(cartCoords, x, coordMaps);
	        % We need to determine the correct coordinate mapping for each
	        % node
    	    coordMaps = computeCoordinateMaps(cartCoords, warpCartCoords);
	        x = x(:);

	        % fminunc only works with data types of class double
    	    [derivWarp, f, exitflag, output] = fminunc(@computeObjective, cast(derivWarp, 'double'), options);
        	displayLogItem(['Results of optimization at resolution #', num2str(res), ':'], logFile);
            	displayLogItem(sprintf('\t%s', ['Number of iterations:  ', num2str(output.iterations)]), logFile);            
	            displayLogItem(sprintf('\t%s', ['Number of function evaluations:  ', num2str(output.funcCount)]), logFile);
    	        displayLogItem(sprintf('\t%s', ['Algorithm used:  ', output.algorithm]), logFile);
        	    displayLogItem(sprintf('\t%s', ['Exit message:  ', output.message]), logFile);
        	    
        	x(derivWarpInds) = derivWarp;
		else 	       
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
		end
        
        
        % Now we have to trim neighbors from the outermost resolution
        % First, decrease the neighbors count for each voxel by the number
        % of voxel neighbors in the outermost resolution
        totalNbrs = totalNbrs - resNbrSizes(end, :);
        % Next, remove the neighbor count for the outermost resolution
        resNbrSizes = resNbrSizes(1:end-1, :);
        
        updTotalNbrs = updTotalNbrs - updResNbrSizes(end, :);
        updResNbrSizes = updResNbrSizes(1:end-1, :);
        
        % We could also re-trim the neighbors matrix
        maxNbrs = cast(max(max(totalNbrs), max(updTotalNbrs)), 'double');
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

		displayLogItem(sprintf('\t%s', ['Iteration number:  ', num2str(optimVals.iteration)]), logFile);
			displayLogItem(sprintf('\t\t%s', ['Function value:  ', num2str(optimVals.fval)]), logFile);
			displayLogItem(sprintf('\t\t%s', ['Step size:  ', num2str(optimVals.stepsize)]), logFile);
			
	        if hasU
		        corrAvg = 0;
    		    for j = 1:N
        		    corrAvg = corrAvg + W2TU1(j, :)*Q(:, j)/N;
		        end

    		    displayLogItem(sprintf('\t\t%s', ['Average inter-subject correlation:  ', num2str(corrAvg)]), logFile);        
    		end
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

    if strcmp(derivCheck, 'on')       
    	xBig = zeros(2*N, 1);
    	xBig(derivWarpInds) = x;
    	x = xBig;
    	clear xBig;
    	
    	reshape(x(derivWarpInds), 2, numDerivCheckNodes)
    end
    
	displayLogItem(sprintf('\t%s', 'Entering compute objective function...'), logFile);
    x = cast(reshape(x, 2, N), opsCoordsType);
	pause(0.01);

	displayLogItem(sprintf('\t\t%s', 'Computing warp cart coords from spherical warp...'), logFile);
	    % First, compute the updated warped cartesian coordinates
    	warpCartCoords = computeWarpCoordsFromSphericalWarp(cartCoords, x, coordMaps);
	displayLogItem(sprintf('\t\t%s', 'Completed computing warp cart coords from spherical warp'), logFile);

	displayLogItem(sprintf('\t\t%s', 'Computing the coordinate maps for each node...'), logFile);
	    % We need to determine the correct coordinate mapping for each node
	    coordMaps = computeCoordinateMaps(cartCoords, warpCartCoords);
	displayLogItem(sprintf('\t\t%s', 'Completed computing the coordinate maps for each node'), logFile);

	displayLogItem(sprintf('\t\t%s', 'Computing the updated neighbors for each node...'), logFile);
	    % Update the neighbors data structure 
	    [updNbrs, updResNbrSizes, updTotalNbrs] = updateNeighborResolutions(cartCoords, warpCartCoords, nbrs, resNbrSizes, totalNbrs, updNbrs, updResNbrSizes, updTotalNbrs);
	displayLogItem(sprintf('\t\t%s', 'Completed computing the updated neighbors for each node'), logFile);
        
	displayLogItem(sprintf('\t\t%s', 'Computing the spherical coords from cartesian...'), logFile);
	    % Now compute the warped locations in spherical coordinates
    	warpSpherCoords = computeSphericalFromCartesian(warpCartCoords,	coordMaps);
	displayLogItem(sprintf('\t\t%s', 'Completed computing the spherical coords from cartesian'), logFile);
    
	% Re-initialize all data structures
	Q = zeros(M1, N, opsDataType);
	if nargout >= 2
        % Then we must compute gradient -- allocate space for it
        dAD_dphi = zeros(maxNbrs, N, opsDataType);
        dAD_dtheta = zeros(maxNbrs, N, opsDataType);
        
        g = zeros(2*N, 1, opsCoordsType);
    end
	if nargout >= 3
        % Then we must compute Hessian -- allocate space for it
	    if strcmp(derivCheck, 'on')       
	    	H = sparse(2*numDerivCheckNodes, 2*numDerivCheckNodes);
	    	H = spdiags(ones(2*numDerivCheckNodes, 1), 0, H);
		else 	       
	        H = sparse(2*N, 2*N);
	        H = spdiags(ones(2*N, 1), 0, H);
    	end
    end
    
	if strcmp(compMethod, 'fast')
		% Compute A matrix
		% If nargout >= 2, then also compute dA_dx and dA_dy
		displayLogItem(sprintf('\t\t%s', 'Computing interpolation matrix A...'), logFile);

    	    % Find the locations where no neighbor exists (i.e., updNbrs == -99)
        	neg99Locs = find(updNbrs == -99);
	        fullUpdNbrs = updNbrs;
    	    fullUpdNbrs(neg99Locs) = floor((neg99Locs-1)/maxNbrs);
	        fullUpdNbrsV1S = fullUpdNbrs;
	        fullUpdNbrsV1S(neg99Locs) = N; % (node #N corresponds to row (N+1) in V1S0)

			displayLogItem(sprintf('\t\t\t%s', 'Reshaping cart coord matrices...'), logFile);
	        	% warp cartcoords
	    	    wcc = reshape(repmat(warpCartCoords, maxNbrs, 1), 3, N*maxNbrs);
    	    	% neighbor cartcoords (NON WARPED!!!)
	    	    nc = cartCoords(:, fullUpdNbrs(:)+1);
			displayLogItem(sprintf('\t\t\t%s', 'Completed reshaping cart coord matrices'), logFile);

			displayLogItem(sprintf('\t\t\t%s', 'Computing geodesic distances...'), logFile);
				% Compute geodesic distances between coordinates
    	    	gds = computeGeodesicDistances(wcc, nc);
			displayLogItem(sprintf('\t\t\t%s', 'Completed computing geodesic distances'), logFile);
    	    
			displayLogItem(sprintf('\t\t\t%s', 'Computing interpolation coefficients...'), logFile); 	   	    
	    	    % Avals will be N*maxNbrs x 1
    	    	[AVals, nonZeroLocs] = interp_f(wcc, nc, res, gds, opsDataType, neg99Locs);
			displayLogItem(sprintf('\t\t\t%s', 'Completed computing interpolation coefficients'), logFile); 	   	    
        
			displayLogItem(sprintf('\t\t\t%s', 'Computing Q, and normalizing its columns...'), logFile); 	   	    
    	    A = reshape(AVals, maxNbrs, N);        
    	    clear AVals;
        	ACell = mat2cell(A', ones(1, N));
	        D = zeros(1, N, opsDataType);
    	    % Q(:, j) = V1ST(:, NI)*A;
	        if nargout >= 2
    	        Q_V1S = zeros(N, maxNbrs, opsDataType);
        	end
	        NSplitSize = ceil(N/maxNbrs);
    	    nStart = NSplitSize*(0:(maxNbrs-1));
        	nStart(nStart > N) = [];
			for ns = nStart
    	        ne = min(ns + NSplitSize, N);
        	    NC = (maxNbrs*ns:(maxNbrs*ne - 1)) + 1;
				V1SCell = mat2cell(V1S0(fullUpdNbrsV1S(NC)+1, 1:M1), maxNbrs*ones(1, (ne-ns)));
    	        Q(1:M1, (ns+1):ne) = transpose(cellMultPre(ACell((ns+1):ne), V1SCell));
        	    mags = sqrt(sum(Q(:, (ns+1):ne).^2));
	            inds = (ns+1):ne;
    	        mLocs = find(mags > 1e-8);
        	    dInds = inds(mLocs);
	            D(dInds) = 1./mags(mLocs);
    	        Q(1:M1, (ns+1):ne) = diagMultPost(Q(:, (ns+1):ne), D((ns+1):ne));
            
        	    if nargout >= 2
            	    % Then we should compute Q_V1S as well...
	              QCell = mat2cell(Q(:, (ns+1):ne), M1, ones(1, (ne-ns)));
	              % Q_V1S = Q(:, j)'*V1ST(:, NI)
    	          Q_V1S((ns+1):ne, :) = transpose(cellMultPost(V1SCell', QCell));
	            end    
    	    end
    	    clear ACell QCell;
			displayLogItem(sprintf('\t\t\t%s', 'Completed computing Q, and normalizing its columns'), logFile); 	   	    
        
    	    % Now work on the gradient if necessary
        	if nargout >= 2
            	% Now we must compute gradient information
	            % This computes derivative of interpolation coefficients with
    	        % respect to phi and theta
        	    cM = cast(reshape(repmat(coordMaps, maxNbrs, 1), 1, maxNbrs*N), 'double');

				displayLogItem(sprintf('\t\t\t%s', 'Computing spherCoord1...'), logFile); 
					spherCoord1 = reshape(repmat(warpSpherCoords, maxNbrs, 1), 2, N*maxNbrs);
				displayLogItem(sprintf('\t\t\t%s', 'Completed computing spherCoord1...'), logFile); 	   	    

				displayLogItem(sprintf('\t\t\t%s', 'Computing spherCoord2...'), logFile); 	   	    
					spherCoord2 = zeros(2, N*maxNbrs, opsCoordsType);
					colsArr = cast(fullUpdNbrs(:)', 'double');
					sArr = zeros(2, N*maxNbrs, opsCoordsType);
					sArr(1, :) = 6*colsArr + 2*(cM-1)+1;
					sArr(2, :) = sArr(1, :) + 1;
					spherCoord2(1:2, 1:N*maxNbrs) = sC(sArr);
				displayLogItem(sprintf('\t\t\t%s', 'Completed computing spherCoord2...'), logFile); 	   	    
            
				displayLogItem(sprintf('\t\t\t%s', 'Computing gradient information...'), logFile); 	   	    
        	    % Now compute dA_dphi and dA_dtheta
            	% These are derivatives of interpolation coefficients with
	            % respect to phi and theta
    	        dA_dphi = reshape(df_dphi(spherCoord1, spherCoord2, res, gds), maxNbrs, N);
        	    dA_dphi(neg99Locs) = 0;
            	dA_dtheta = reshape(df_dtheta(spherCoord1, spherCoord2, res, gds), maxNbrs, N);
	            dA_dtheta(neg99Locs) = 0;
    	        
        	    dD_dphi = (D.^2).*sum(Q_V1S'.*dA_dphi);
            	dD_dtheta = (D.^2).*sum(Q_V1S'.*dA_dtheta);
            	clear Q_V1S;
            
	            dAD_dphi = diagMultPost(dA_dphi, D) - diagMultPost(A, dD_dphi);
    	        dAD_dtheta = diagMultPost(dA_dtheta, D) - diagMultPost(A, dD_dtheta);
				displayLogItem(sprintf('\t\t\t%s', 'Completed computing gradient information'), logFile); 	   	    
			end
    elseif strcmp(compMethod, 'slow')    
		% OLD METHOD FOR INTERPOLATING
		% For each j, compute the jth row of A
	    displayLogItem(sprintf('\t\t%s', 'Computing interpolation matrix A!'), logFile);            
		for j = 1:N
			currCoordMap = coordMaps(:, j);
			spherCoords = sC((2*currCoordMap-1):2*currCoordMap, :);%spherCoordsCellArr{currCoordMap};
            
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

%			Q(:, j) = V1ST(:, NI)*A;
			Q(:, j) = transpose(A'*V1S(NI, :));

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

%				Q_V1S = Q(:, j)'*V1ST(:, NI);
				Q_V1S = transpose(V1S(NI, :)*Q(:, j));

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
		% END:  OLD METHOD FOR INTERPOLATING
	else
        error('unknown comp mode');
    end
    
    f = 0;        
    if lambda_intrasubj > 0
    	displayLogItem(sprintf('\t\t%s', 'Computing intra-subject correlation portion of the objective...'), logFile);
    		displayLogItem(sprintf('\t\t\t%s', 'Computing B1S...'), logFile);
	    		B1S = Q*V2;
    		displayLogItem(sprintf('\t\t\t%s', 'Completed computing B1S!'), logFile);
    		
    		displayLogItem(sprintf('\t\t\t%s', 'Computing null-space basis for V2...'), logFile);
    			V2_0 = Q' - V2*B1S';
    		displayLogItem(sprintf('\t\t\t%s', 'Completed computing null-space basis for V2!'), logFile);
    		
    		displayLogItem(sprintf('\t\t\t%s', 'Computing QR factorization of null-space basis for V2...'), logFile);
    			V2_0 = orth_qr(V2_0);
    		displayLogItem(sprintf('\t\t\t%s', 'Completed computing QR factorization of null-space basis for V2!'), logFile);
    		
    		displayLogItem(sprintf('\t\t\t%s', 'Checking if vectors are orthogonal...'), logFile);
    			maxCrossTerm = max(max(abs(V2'*V2_0)));
    			while maxCrossTerm > 3e-5 % changed from 1e-6 for speed
    				displayLogItem(sprintf('\t\t\t\t%s', ['Vectors not orthogonal -- ', num2str(maxCrossTerm)]), logFile);
    				V2_0 = V2_0 - V2*(V2'*V2_0);
    				V2_0 = orth(V2_0);      % changed from orth_qr() to avoid a bug
    				maxCrossTerm = max(max(abs(V2'*V2_0)));
    			end
    		displayLogItem(sprintf('\t\t\t%s', 'Completed checking if vectors are orthogonal!'), logFile);
    		
    		displayLogItem(sprintf('\t\t\t%s', 'Computing B2S...'), logFile);
    			B2S = Q*V2_0;
    		displayLogItem(sprintf('\t\t\t%s', 'Completed computing B2S!'), logFile);
    		
    		displayLogItem(sprintf('\t\t\t%s', 'Computing intra-subject correlation norm...'), logFile);
    			f = f + lambda_intrasubj*(norm(B1S'*B1S - diag(s2_squared), 'fro')^2 + norm(B2S'*B2S, 'fro')^2 + 2*norm(B1S'*B2S, 'fro')^2);
    		displayLogItem(sprintf('\t\t\t%s', 'Completed computing intra-subject correlation norm!'), logFile);
    		
    	displayLogItem(sprintf('\t\t%s', 'Completed computing intra-subject correlation portion of the objective!'), logFile);
    end
    
    if strcmp(regularization.mode, 'metric_and_areal_distortion')
        % We need to compute the areal and metric distortion
        if regularization.lambda_metric > 0
            displayLogItem(sprintf('\t\t%s', 'Computing metric distortion penalty...'), logFile);
            
            if nargout >= 2
                [mdDiff, dmdDiffs_dphi, dmdDiffs_dtheta] = computeMetricTerms(origNbrs, origNumNbrs, warpCartCoords, coordMaps, regularization.metricDistances, res);
                g(1:2:end) = g(1:2:end) + 4*regularization.lambda_metric*dmdDiffs_dphi;
                g(2:2:end) = g(2:2:end) + 4*regularization.lambda_metric*dmdDiffs_dtheta;
            else
                mdDiff = computeMetricTerms(origNbrs, origNumNbrs, warpCartCoords, coordMaps, regularization.metricDistances, res);
            end
            
            f = f + regularization.lambda_metric*mdDiff;
            displayLogItem(sprintf('\t\t%s', 'Completed computing metric distortion penalty!'), logFile);
        end
        
        if regularization.lambda_areal > 0
        	displayLogItem(sprintf('\t\t%s', 'Computing areal distortion penalty...'), logFile);
            % Now compute the areal term
            if nargout >= 2
                [triArea, dareal_dphi, dareal_dtheta] = computeArealTerms(regularization.triangles, warpCartCoords, coordMaps, regularization.orientedAreas, regularization.orientedNormals);
                g(1:2:end) = g(1:2:end) + 2*regularization.lambda_areal*dareal_dphi;
                g(2:2:end) = g(2:2:end) + 2*regularization.lambda_areal*dareal_dtheta;
            else
                triArea = computeArealTerms(regularization.triangles, warpCartCoords, coordMaps, regularization.orientedAreas, regularization.orientedNormals);                
            end
            f = f + regularization.lambda_areal*triArea;
            displayLogItem(sprintf('\t\t%s', 'Completed computing areal distortion penalty!'), logFile);
        end
    end


    if nargout >= 2
		if strcmp(compMethod, 'fast')
	    	displayLogItem(sprintf('\t\t%s', 'Computing gradient...'), logFile);
    	    if lambda_intrasubj > 0
        	    Mpost = 4*((B1S*B1S' + B2S*B2S')*Q - B1S*V2S2T);
            	MpostCell = mat2cell(Mpost, M1, ones(1, N));
	            dS_dAD = zeros(maxNbrs, N, opsDataType);
    	        for ns = nStart
        	        ne = min(ns + NSplitSize, N);
            	    NC = (maxNbrs*ns:(maxNbrs*ne - 1)) + 1;
                	V1SCell = mat2cell(V1S0(fullUpdNbrsV1S(NC)+1, 1:M1), maxNbrs*ones(1, (ne-ns)));
	                dS_dAD(1:maxNbrs, (ns+1):ne) = cellMultPost(V1SCell', MpostCell((ns+1):ne));
    	        end
    	        clear V1SCell;
        	end
	        g(1:2:end) = g(1:2:end) + transpose(lambda_intrasubj*sum(dS_dAD.*dAD_dphi));
    	    g(2:2:end) = g(2:2:end) + transpose(lambda_intrasubj*sum(dS_dAD.*dAD_dtheta));
        	displayLogItem(sprintf('\t\t%s', 'Completed computing gradient!'), logFile);
		elseif strcmp(compMethod, 'slow')
			% OLD METHOD FOR COMPUTING GRADIENT!!!
			if lambda_intrasubj > 0
				Mpre1 = V1S*(B1S*B1S' + B2S*B2S');
				Mpre2 = V1S*B1S;
            end

            count = 1;
            for p = 1:N
                locLength = locsLength(p);
                locs = Alocs(1:locLength, p);

                if lambda_intrasubj > 0
                    dS_dAD = 4*(Mpre1(locs, :)*Q(:, p) - Mpre2(locs, :)*V2S2T(:, p));
                    g(count) = g(count) + lambda_intrasubj*dS_dAD'*dAD_dphi(1:locLength, p);
                    g(count + 1) = g(count + 1) + lambda_intrasubj*dS_dAD'*dAD_dtheta(1:locLength, p);
                end

                if lambda_intersubj > 0
                    dS_dAD = W2TU1(p, :)*V1S(locs, :)';
                    g(count) = g(count) - lambda_intersubj*dS_dAD*dAD_dphi(1:locLength, p);
                    g(count + 1) = g(count + 1) - lambda_intersubj*dS_dAD*dAD_dtheta(1:locLength, p);
                end

                count = count + 2;
            end
			displayLogItem(sprintf('\t\t%s', 'Completed computing gradient!'), logFile);
			% END OLD METHOD FOR COMPUTING GRADIENT!!!       
		else
			error('unknown comp mode');
		end
	    if strcmp(derivCheck, 'on')  
   			g = g(derivWarpInds);
    	end
    end
    
    if strcmp(derivCheck, 'on')
    	displayLogItem(sprintf('\t\t%s', ['Objective f:  ', num2str(f)]), logFile);
    end

	displayLogItem(sprintf('\t%s', 'Exiting compute objective function!'), logFile);   
    % Force f to be double and make sure f and g are real
    % to match requirements of optimization functions in Matlab
    if isa(f,'double') == 0
        f = double(f);
    end
    g = real(g);
    f = real(f);
end
%% ** END FUNCTION computeObjective
end
%% ** END connnorm_register
