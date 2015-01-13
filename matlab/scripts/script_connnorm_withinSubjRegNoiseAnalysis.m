function warps = script_connnorm_withinSubjRegNoiseAnalysis(subjects, hems, experiment, lambda, noise_sigma)
% FUNCTION warps = script_connnorm_withinSubjReg(subjects, hems, experiment, lambda, noise_sigma)
% *** INPUT ARGUMENTS ***
%   subjects:  a cell array of strings specifying the subject names
%       Ex:  {'cb', 'dm', 'hj', 'kd', 'kl', 'mh', 'ph', 'rb', 'se', 'sm'}
%   hems:  a cell array of strings specifying the hemispheres to align
%       Ex:  {'lh', 'rh'}
%   experiment:  name of the functional experiment (e.g., 'movie5min_P1')
%   mmRes:  the resolution of the cortical mesh, in mm (i.e., 2)
%   fwhm:  the fwhm for SurfSmooth (i.e., 0)
%   fs_surf:  the freesurfer surface file -- (i.e., "sphere", or
%       "sphere.reg")
%   outDir:  the specified output directory 
%   lambda_areal:  parameter weighting for the areal distortion penalty
%       term
%   lambda_metric:  parameter weighting for the metric distortion penalty
%       term
%   maxResolution:  multiplier specifying the maximum resolution for multi-resolution alignment
%       The alignment begins at a resolution of maxResolution*mmRes, where
%       mmRes is specified above as the initial cortical mesh resolution
%   logFile (optional):  filename to store log items (relative to outDir)
%
% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory


	mmRes = 2;
	fwhm = 0;
	fs_surf = 'sphere.reg';
	outDir = '../surf_data/connnorm_out_data/';
	maxResolution = 3;
	logFile = 'connnorm_log.log';
	pctEnergy = 95;
	
	lambda_areal = lambda;
	lambda_metric = lambda;
	
    saveDataBytes = 4; % Number of bytes of datatype to store data (e.g. float32 = 4 bytes)
    saveDataType = 'single';
    opsDataBytes = 8; % Number of bytes to perform all operations in (e.g. single precision = 4 bytes)
    opsDataType = 'single';
    opsCoordsType = 'double'; % Format for coords structures
    saveCoordsType = 'single';
    nsubj = length(subjects); % Number of subjects
    nhem = length(hems); % Number of hemispheres
    
    displayLogItem('Beginning conn-norm algorithm...', logFile);

%% Setup output directory structure
    if outDir(end) ~= '/'
        outDir = [outDir, '/'];
    end
    
    % The algorithm creates a subdirectory under outDir to save all files related to this alignment
    alignDir = [outDir, 'alignment_withinSubjRegNoiseAnalysisPreNormalized_', experiment, '_noise_', num2str(noise_sigma), '_', num2str(lambda_areal), '/'];
    if isdir(alignDir) == 0
        mkdir(alignDir);
    else
    	error('directory already exists!');
    end
    
    % logDir is where all log files are saved
    logDir = [alignDir, 'log/'];
        if isdir(logDir) == 0
            mkdir(logDir);
        end
    % warpsDir is where all warps are saved
    warpsDir = [alignDir, 'warps/'];
        if isdir(warpsDir) == 0
            mkdir(warpsDir);
        end
        
    if strlen(logFile)
    	logFile = [logDir, logFile];
   	else
    	logFile = '';	
   	end
    
    determineStartPoint = 0;

    tmpDir = [alignDir, 'tmp/'];
    if isdir(tmpDir) == 0
        mkdir(tmpDir);
    end

    for hemNum = 1:nhem
        hem = hems{hemNum};

        hemOutDirs{hemNum} = [warpsDir, hem, '/'];
        if isdir(hemOutDirs{hemNum}) == 0
            mkdir(hemOutDirs{hemNum});
        else
            % It appears that this run has been started in the past.  Set flag
            % to determine starting point
            determineStartPoint = 1;
        end
    end
%% END Setup output directory structure


    displayLogItem('Beginning connnorm algorithm...', logFile); 
    
    displayLogItem('Miscellaneous parameters:', logFile);
	    displayLogItem(sprintf('\t%s', ['Experiment:  ', experiment]), logFile);
    	displayLogItem(sprintf('\t%s', ['MM resolution:  ', num2str(mmRes)]), logFile);
    	displayLogItem(sprintf('\t%s', ['FWHM:  ', num2str(fwhm)]), logFile);
    	displayLogItem(sprintf('\t%s', ['lambda_areal:  ', num2str(lambda_areal)]), logFile);
    	displayLogItem(sprintf('\t%s', ['lambda_metric:  ', num2str(lambda_metric)]), logFile);
    	displayLogItem(sprintf('\t%s', ['Max resolution:  ', num2str(maxResolution)]), logFile);
	

%% Parse the surface file and setup coordinates structure
    displayLogItem(['Parsing the surface file from directory ', outDir, '...'], logFile);
    surfFile = [outDir, 'standard', num2str(mmRes), 'mm_', fs_surf, '.asc'];
    [numNodesPerHem, coords] = parseSurfaceFile(surfFile, opsCoordsType, saveCoordsType);
    numNodes = numNodesPerHem*nhem;
    numTrisPerHem = size(coords.triangles, 2);
    numTris = nhem*numTrisPerHem;
    % Normalize the coordinates to the unit sphere
    coords.cartCoords = normalizeCartesianCoordinates(coords.cartCoords);
    if nhem > 1
	    % Replicate everything for all hemispheres
        coords = generateMultiHemisphereCoordsDataset(coords, nhem);
    end
    displayLogItem('Completed parsing the surface file', logFile);

	noise_res = 3;
	if noise_sigma > 0
	    [nbrs, resNbrSizes, totalNbrs] = computeNeighborResolutions(coords.neighbors, noise_res);
	end
%% END Parse the surface file and setup coordinates structure
	

    % Determine the number of time-points
    tLen = getNumTimepoints(getOrigTimeSeriesFilename(subjects{1}, hems{1}));
	displayLogItem(['Number of nodes per hemisphere:  ', num2str(numNodesPerHem)], logFile);
	displayLogItem(['Number of triangles per hemisphere:  ', num2str(numTrisPerHem)], logFile);
	displayLogItem(['Number of timepoints:  ', num2str(tLen)], logFile);
	displayLogItem(['Number of hemispheres:  ', num2str(nhem)], logFile);
	
%% Setup regularization
    regularization.mode = 'metric_and_areal_distortion';
    displayLogItem('Using metric and areal distortion as warp regularization', logFile);
    displayLogItem(sprintf('\tLambda areal:  %f', lambda_areal), logFile);
    displayLogItem(sprintf('\tLambda metric:  %f', lambda_metric), logFile);
    regularization.lambda_areal = numNodes*lambda_areal;
    regularization.lambda_metric = numNodes*lambda_metric;
%% END Setup regularization


for subjNum = 1:nsubj
	subj = subjects{subjNum};
	
	T1 = loadOrigTimeseries(subj);
	T1 = normalizeDataset(T1);
	if noise_sigma > 0
		% Generate additive noise
		e = noise_sigma*rand(size(T1));
		
		% Temporally smooth
		H = fspecial('gaussian', [5, 1]);
		e = filter(H, 1, e);
		
		% Spatially smooth
		[eu,es,ev] = svd(e,'econ');
		es=diag(es);
		clear e;
		[ev,es,eu] = blurDataset(ev, es, eu, coords.cartCoords, nbrs, totalNbrs, noise_res);
		e = eu*diag(es)*ev';
		clear eu es ev;
		
		T2 = T1 + e;
		clear e;
	else
		T2 = T1;
	end
	T2 = normalizeDataset(T2);
		
	[u,s,v] = svd(T1, 'econ');
	s = diag(s);
	clear u;
	p = findPercentile(s, pctEnergy);
	s = s(1:p);
	v = v(:, 1:p);
	displayLogItem(sprintf('\t%s', ['Keeping ', num2str(p), ' components']), logFile);
	sub1stHalf.sigma = s;
	sub1stHalf.V = v;
	clear T1 s v;
	
	[u,s,v] = svd(T2, 'econ');
	s = diag(s);
	clear u;
	p = findPercentile(s, pctEnergy);
	s = s(1:p);
	v = v(:, 1:p);
	displayLogItem(sprintf('\t%s', ['Keeping ', num2str(p), ' components']), logFile);
	sub2ndHalf.sigma = s;
	sub2ndHalf.V = v;
	clear T2 s v;
	
	% Set the output file prefix
	alignLogFile = [logDir, 'log_alignment_noise_subjNum_', num2str(subjNum), '_subj_', subj, '.txt'];
	% Now run the algorithm
	displayLogItem(['Performing alignment noise of subject #', num2str(subjNum), ' (', subj, ')'], logFile);
	displayLogItem(['A log of this alignment will be saved to ', alignLogFile], logFile);
            
    % This aligns the 1st half to the second half
	warp1stTo2nd = connnorm_register(sub1stHalf, sub2ndHalf, coords, zeros(3, numNodes, opsCoordsType), regularization, maxResolution, alignLogFile);
		displayLogItem(['Completed alignment noise of subject #', num2str(subjNum), ' (', subj, ')'], logFile);

	clear sub1stHalf sub2ndHalf;
	% Now we need to output the warp file
	saveWarp(warp1stTo2nd, subj);
	
%	alignLogFile = [logDir, 'log_alignment_2ndHalf_to_1stHalf_subjNum_', num2str(subjNum), '_subj_', subj, '.txt'];
%	displayLogItem(['Performing alignment (2nd half to 1st half) of subject #', num2str(subjNum), ' (', subj, ')'], logFile);
%	displayLogItem(['A log of this alignment will be saved to ', alignLogFile], logFile);
%	warp2ndTo1st = connnorm_register(sub2ndHalf, sub1stHalf, coords, zeros(3, numNodes, opsCoordsType), regularization, maxResolution, alignLogFile);
%		displayLogItem(['Completed alignment (2nd half to 1st half) of subjects #', num2str(subjNum), ' (', subj, ')'], logFile);
		
	% Now we need to output the warp file
%	saveWarp(warp2ndTo1st, subj, '2ndto1st');
end
displayLogItem('Completed alignment', logFile);
    





%% FUNCTION getOrigTimeSeriesFilename
    function fName = getOrigTimeSeriesFilename(subj, hem)
        fName = [outDir, subj, '_', hem, '_', num2str(mmRes), 'mm_fwhm', num2str(fwhm), '_', experiment, '_on_', fs_surf, '.niml.dset'];
    end
%% END FUNCTION getOrigTimeSeriesFilename


%% FUNCTION loadOrigTimeseries
    function T = loadOrigTimeseries(subj)
        T = zeros(tLen, nhem*numNodesPerHem, opsDataType);

        n = 1;
        for hemNum = 1:nhem
            hem = hems{hemNum};
            T(:, n:(n+numNodesPerHem-1)) = read_niml_binary(getOrigTimeSeriesFilename(subj, hem), opsDataType, saveDataType);

            n = n + numNodesPerHem;
        end
    end
%% END FUNCTION loadOrigTimeseries


%% FUNCTION getWarpFilename
    function warpFilename = getWarpFilename(subj, hem)
        warpFilename = [tmpDir, 'warp_', subj, '_', hem, '.bin'];
    end
%% END FUNCTION getWarpFilename


%% FUNCTION saveWarp
    function saveWarp(warp, subj)
        n = 1;
        for hemNum = 1:nhem
            hem = hems{hemNum};

            fp = fopen(getWarpFilename(subj, hem), 'wb');
            fwrite(fp, warp(:, n:(n+numNodesPerHem-1)), saveCoordsType);
            fclose(fp);

            n = n + numNodesPerHem;
        end
    end
%% END FUNCTION saveWarp


end
%% END script_connnorm_withinSubjReg