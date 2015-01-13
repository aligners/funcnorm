function warps = script_connnorm_acrossSubjReg(subjects, hems, varargin)
% FUNCTION warps = script_connnorm_acrossSubjReg(subjects, hems, experiment, mmRes, fwhm, fs_surf, outDir, ...
%   lambda_areal, lambda_metric, maxResolution [, logFile])
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


    if nargin == 10 || nargin == 11
        % FUNCTION warps = connnorm(subjects, hems, experiment, mmRes, fwhm, fs_surf, outDir, ...
        %   lambda_areal, lambda_metric, maxResolution [, logFile])
        
        % Then don't run funcnorm_prepare.sh
        runFuncnormPrepare = 0;
        
        experiment = varargin{1};
        mmRes = varargin{2};
        fwhm = varargin{3};
        fs_surf = varargin{4};
        outDir = varargin{5};
        lambda_areal = varargin{6};
        lambda_metric = varargin{7};
        maxResolution = varargin{8};
        
        if nargin == 10
            logFile = '';
        else
            logFile = varargin{9};
        end
    else
        error('Improper call to script_connnorm_withinSubjReg.m  Type ''help script_connnorm_withinSubjReg'' for details.');
    end
    
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
    alignDir = [outDir, 'alignment_acrossSubjReg_', experiment, '_', num2str(lambda_areal), '/'];
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
    
%    % Check if multithreading is on...
%    mCT = maxNumCompThreads;
%    if mCT == 1
%    	cTRep = input('Multithreading is currently set to off.  It is strongly recommended that it be turned on.  Would you like to turn it on? Y/N [Y]:  ', 's');
%    	if isempty(cTRep)
%    		cTRep = 'Y';
%    	end
%    	
%    	if cTRep == 'Y';
%    		mCT = maxNumCompThreads('automatic');
%		end
%	end
%	displayLogItem(['Number of computational threads:  ', num2str(mCT)], logFile);
    
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

crossSubjMat = zeros(nsubj, nsubj);
for runNum = 1:nsubj
	% Generate a random subject 1 & 2
	while 1
		subjNum1 = round(0.5 + nsubj*rand(1));
		subjNum2 = round(0.5 + nsubj*rand(1));
		if subjNum1 == subjNum2
			continue;
		end
		if crossSubjMat(subjNum1, subjNum2) == 1
			continue;
		end
		break;
	end
	crossSubjMat(subjNum1, subjNum2) = 1;
	crossSubjMat(subjNum2, subjNum1) = 1;

	subj1 = subjects{subjNum1};
	subj2 = subjects{subjNum2};
	
	T1 = loadOrigTimeseries(subj1);
	t1size = ceil(size(T1,1)/2);
	T2 = loadOrigTimeseries(subj2);

	T1 = normalizeDataset(T1(1:t1size,:));
	T2 = normalizeDataset(T2((t1size+1):end,:));
	displayLogItem(['Number of timepoints in 1st half:  ', num2str(t1size)], logFile);
	displayLogItem(['Number of timepoints in 2nd half:  ', num2str(size(T2, 1))], logFile);
	
	[u,s,v] = svd(T1, 'econ');
	s = diag(s);
	clear u;
	sub1stHalf.sigma = s;
	sub1stHalf.V = v;
	clear u s v;
	[u,s,v] = svd(T2, 'econ');
	s = diag(s);
	sub2ndHalf.sigma = s;
	sub2ndHalf.V = v;
	clear T u s v;
	
	% Set the output file prefix
	alignLogFile = [logDir, 'log_alignment_subj_', subj1, '_1stHalf_to_subj_', subj2, '_2ndHalf.txt'];
	% Now run the algorithm
	displayLogItem(['Performing alignment (1st half to 2nd half) of subject #', num2str(subjNum1), ' to subject #', num2str(subjNum2)], logFile);
	displayLogItem(['A log of this alignment will be saved to ', alignLogFile], logFile);
            
    % This aligns the 1st half to the second half
	warp1stTo2nd = connnorm_register(sub1stHalf, sub2ndHalf, coords, zeros(3, numNodes, opsCoordsType), regularization, maxResolution, alignLogFile);
		displayLogItem(['Completed alignment'], logFile);

	% Now we need to output the warp file
	saveWarp(warp1stTo2nd, runNum);
	
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
    function warpFilename = getWarpFilename(runNum, hem)
        warpFilename = [tmpDir, 'warp_runNum_', num2str(runNum), '_', hem, '.bin'];
    end
%% END FUNCTION getWarpFilename


%% FUNCTION saveWarp
    function saveWarp(warp, runNum)
        n = 1;
        for hemNum = 1:nhem
            hem = hems{hemNum};

            fp = fopen(getWarpFilename(runNum, hem), 'wb');
            fwrite(fp, warp(:, n:(n+numNodesPerHem-1)), saveCoordsType);
            fclose(fp);

            n = n + numNodesPerHem;
        end
    end
%% END FUNCTION saveWarp


end
%% END script_connnorm_withinSubjReg