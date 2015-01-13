function warps = script_connnorm_localMinAnalysis(subjects, hems, varargin)
% FUNCTION warps = script_connnorm_localMinAnalysis(subjects, hems, experiment, mmRes, fwhm, fs_surf, outDir, ...
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
    
    inpWarpsDir = [outDir, 'alignment_', experiment, '/tmp/'];
    
    % The algorithm creates a subdirectory under outDir to save all files related to this alignment
    alignDir = [outDir, 'alignment_localMinAnalysis_', experiment, '_', num2str(lambda_areal), '_', num2str(maxResolution), '/'];
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
        
    if length(logFile)
    	logFile = [logDir, logFile];
   	else
    	logFile = '';	
   	end
    
    tmpDir = [alignDir, 'tmp/'];
    if isdir(tmpDir) == 0
        mkdir(tmpDir);
    end

    for hemNum = 1:nhem
        hem = hems{hemNum};

        hemOutDirs{hemNum} = [warpsDir, hem, '/'];
        if isdir(hemOutDirs{hemNum}) == 0
            mkdir(hemOutDirs{hemNum});
        end
    end
%% END Setup output directory structure


    displayLogItem('Beginning connnorm algorithm...', logFile); 
    
    % Check if multithreading is on...
    mCT = maxNumCompThreads;
    if mCT == 1
    	cTRep = input('Multithreading is currently set to off.  It is strongly recommended that it be turned on.  Would you like to turn it on? Y/N [Y]:  ', 's');
    	if isempty(cTRep)
    		cTRep = 'Y';
    	end
    	
    	if cTRep == 'Y';
    		mCT = maxNumCompThreads('automatic');
		end
	end
	displayLogItem(['Number of computational threads:  ', num2str(mCT)], logFile);
    
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


for subjNum = 1:nsubj
	subj = subjects{subjNum};
	
	T = centerDataset(loadOrigTimeseries(subj));
	[u,s,v] = svd(T, 'econ');
	s = diag(s);
	clear u;
	sub.sigma = s;
	sub.V = v;
	clear u s v;
	
    fp = fopen([inpWarpsDir, 'warp_', subj, '_lh_0.bin'], 'rb');
    warp(1:3, 1:numNodesPerHem) = fread(fp, [3, numNodesPerHem], 'float32');fclose(fp);
    fp = fopen([inpWarpsDir, 'warp_', subj, '_rh_0.bin'], 'rb');
    warp(1:3, (numNodesPerHem+1):numNodes) = fread(fp, [3, numNodesPerHem], 'float32');fclose(fp);
    
	% Set the output file prefix
	alignLogFile = [logDir, 'log_alignment_subjNum_', num2str(subjNum), '_subj_', subj, '.txt'];
	% Now run the algorithm
	displayLogItem(['Performing alignment of subject #', num2str(subjNum), ' (', subj, ')'], logFile);
	displayLogItem(['A log of this alignment will be saved to ', alignLogFile], logFile);
	warp = connnorm_register(sub, sub, coords, warp, regularization, maxResolution, alignLogFile);
    displayLogItem(['Completed alignment of subject #', num2str(subjNum), ' (', subj, ')'], logFile);

	% Now we need to output the warp file
	saveWarp(warp, subj);
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