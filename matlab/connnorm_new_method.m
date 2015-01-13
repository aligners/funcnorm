function warps = connnorm(subjects, hems, varargin)
% FUNCTION warps = connnorm(subjects, hems, varargin)
% *** INPUT ARGUMENTS ***
%   subjects:  a cell array of strings specifying the subject names
%       Ex:  {'cb', 'dm', 'hj', 'kd', 'kl', 'mh', 'ph', 'rb', 'se', 'sm'}
%   hems:  a cell array of strings specifying the hemispheres to align
%       Ex:  {'lh', 'rh'}
%
%   Additional input arguments depend on whether the data should first be
%   prepared with the funcnorm_prepare.sh script
%
%
% FUNCTION warps = connnorm(subjects, hems, experiment, mmRes, fwhm, fs_surf, outDir, ...
%   lambda_areal, lambda_metric, pctEnergy, maxResolution, numPasses [, logFile [, tLen ]])
% Use this form only if the data has already been prepared using the
% funcnorm_prepare.sh script.  
% If funcnorm_prepare.sh has not yet been run on the data, see the function
% prototype below
% *** INPUT ARGUMENTS ***
%   subjects:  see above
%   hems:  see above
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
%   pctEnergy:  the percentage of fMRI time-series energy to preserve in the PCA
%       dimensionality reduction step
%   maxResolution:  multiplier specifying the maximum resolution for multi-resolution alignment
%       The alignment begins at a resolution of maxResolution*mmRes, where
%       mmRes is specified above as the initial cortical mesh resolution
%   numPasses:  the number of alignment passes through the subjects
%   logFile (optional):  filename to store log items (relative to outDir)
%
%
% FUNCTION warps = connnorm(subjects, hems, experiment, volDataDir, ...
%   suffix, extensions, inpSubjDir, saveSubjDir, preprocess_movie, ...
%   mmRes, fwhm, fs_surf, outDir, lambda_areal, lambda_metric, maxResolution, ... 
%   numPasses [, logFile])
% Use this form if the data has not yet been prepared with the
% funcnorm_prepare.sh script
% *** INPUT ARGUMENTS ***
%   subjects:  see above
%   hems:  see above
%   experiment:  name of the functional experiment (e.g., 'movie5min_P1')
%   volDataDir:  directory location to the volume data
%       Ex:  '/data/subjects_in/sample_data/movie/'
%   suffix:  suffix to volume names
%       Ex:  '+orig';
%   extensions:  the extension to the anatomical and functional volumes
%       Ex:  'BRIK.gz'
%   inpSubjDir:  input subjects directory
%       Ex:  '/data/subjects_in/sample_data/SUBJECTS_DIR'
%   saveSubjDir:  directory to save generated surfaces and intermediate
%       results
%   preprocess_movie:  a flag (0 or 1) specifying whether to perform movie
%       preprocessing on the functional data
%   mmRes:  the resolution of the cortical mesh, in mm (i.e., 2)
%   fwhm:  the fwhm for SurfSmooth (i.e., 0)
%   fs_surf:  the freesurfer surface file -- (i.e., "sphere", or
%       "sphere.reg")
%   outDir:  the specified output directory 
%   lambda_areal:  parameter weighting for the areal distortion penalty
%       term
%   lambda_metric:  parameter weighting for the metric distortion penalty
%       term
%   pctEnergy:  the percentage of time-series energy to preserve in the PCA
%       dimensionality reduction step
%   maxResolution:  multiplier specifying the maximum resolution for multi-resolution alignment
%       The alignment begins at a resolution of maxResolution*mmRes, where
%       mmRes is specified above as the initial cortical mesh resolution
%   numPasses:  the number of alignment passes through the subjects
%   logFile (optional):  filename to store log items (relative to outDir)
%
% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory


    if nargin == 12 || nargin == 13 || nargin == 14
        % FUNCTION warps = connnorm(subjects, hems, experiment, mmRes, fwhm, fs_surf, outDir, ...
        %   lambda_areal, lambda_metric, pctEnergy, maxResolution, numPasses [, logFile])
        
        % Then don't run funcnorm_prepare.sh
        runFuncnormPrepare = 0;
        
        experiment = varargin{1};
        mmRes = varargin{2};
        fwhm = varargin{3};
        fs_surf = varargin{4};
        outDir = varargin{5};
        lambda_areal = varargin{6};
        lambda_metric = varargin{7};
        pctEnergy = varargin{8};
        maxResolution = varargin{9};
        numPasses = varargin{10};
        
        if nargin == 12
        	tLen = -1;
            logFile = '';
        elseif nargin == 13
			tLen = -1;
        	logFile = varargin{11};
        else
            logFile = varargin{11};
        	tLen = varargin{12};
        end
    elseif nargin == 18 || nargin == 19
        % FUNCTION warps = connnorm(subjects, hems, experiment, volDataDir, ...
        %   suffix, extensions, inpSubjDir, saveSubjDir, preprocess_movie, ...
        %   mmRes, fwhm, fs_surf, outDir, lambda_areal, lambda_metric, pctEnergy, ...
        %   maxResolution, numPasses [, logFile])
        
        % Then run funcnorm_prepare.sh first
        runFuncnormPrepare = 1;
        
        experiment = varargin{1};
        volDataDir = varargin{2};
        if volDataDir(end) ~= '/'
            volDataDir = [volDataDir, '/'];
        end
        suffix = varargin{3};
        extensions = varargin{4};
        inpSubjDir = varargin{5};
        if inpSubjDir(end) == '/'
            inpSubjDir = inpSubjDir(1:(end-1));
        end
        saveSubjDir = varargin{6};
        if saveSubjDir(end) ~= '/'
            saveSubjDir = [saveSubjDir, '/'];
        end
        preprocess_movie = varargin{7};
        mmRes = varargin{8};
        fwhm = varargin{9};
        fs_surf = varargin{10};
        outDir = varargin{11};
        lambda_areal = varargin{12};
        lambda_metric = varargin{13};
        pctEnergy = varargin{14};
        maxResolution = varargin{15};
        numPasses = varargin{16};
        
        if nargin == 18
        	logFile = '';
        else
            logFile = varargin{17};
        end
    else
        error('Improper call to connnorm.m  Type ''help connnorm'' for details.');
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
    if tLen ~= -1
	    alignDir = [outDir, 'alignment_', experiment, '_', num2str(lambda_areal), '_', num2str(pctEnergy), '_tLen', num2str(tLen), '/'];
	else
		alignDir = [outDir, 'alignment_', experiment, '_', num2str(lambda_areal), '_', num2str(pctEnergy), '/'];
	end
	
    if isdir(alignDir) == 0
        mkdir(alignDir);
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
    	displayLogItem(sprintf('\t%s', ['Pct Energy:  ', num2str(pctEnergy)]), logFile);
    	displayLogItem(sprintf('\t%s', ['Max resolution:  ', num2str(maxResolution)]), logFile);
    	displayLogItem(sprintf('\t%s', ['Num passes:  ', num2str(numPasses)]), logFile);
	

%% Prepare the data if necessary
    if runFuncnormPrepare == 1
        funcnorm_prepare_wrapper(subjects, hems, experiment, volDataDir, suffix, ...
            extensions, inpSubjDir, saveSubjDir, preprocess_movie, mmRes, fwhm, ...
            fs_surf, outDir, logFile)
    else
    	displayLogItem('Not running funcnorm_prepare script', logFile);
    end
%% END Prepare the data if necessary
    

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
    if tLen == -1
	    tLen = getNumTimepoints(getOrigTimeSeriesFilename(subjects{1}, hems{1}));
	end
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

	atlasProjMethod = 'keepSpecified'; % Can set to 
	absMaxAtlasComponents = 600;
	if strcmp(atlasProjMethod, 'keepAll')
		displayLogItem('There will not be any projection of the atlas principal components.  This will cause the algorithm to be slower, but more accurate', logFile);
	elseif strcmp(atlasProjMethod, 'keepSameAsSubject')
		displayLogItem('The atlas principal components will be projected.  This will cause the algorithm to be faster, but less accurate', logFile);
	elseif strcmp(atlasProjMethod, 'keepSpecified')
		displayLogItem(['The registration will keep ', num2str(absMaxAtlasComponents), ' components in the atlas'], logFile);
	else
		displayLogItem(['ERROR:  Unknown atlasProjMethod:  ', atlasProjMethod], logFile);
		error(['ERROR:  Unknown atlasProjMethod:  ', atlasProjMethod]);
	end

%% Data preprocessing
    displayLogItem('Pre-processing the datasets...', logFile);      
    for subjNum = 1:nsubj
        subj = subjects{subjNum};
        
        [V_fName, s_fName] = getOrigProjectedDataFilenames(subj);
        if (exist(V_fName)==2) && (exist(s_fName)==2)
        	displayLogItem(['SVD of data already found, delete to regenerate'], logFile);
			continue;
        end
        
        % First, load the data and center it
        T = loadOrigTimeseries(subj);
        T = centerDataset(T);
        
        % Compute the singular vectors/values of the time-series matrix
		displayLogItem(sprintf('\t%s', ['Computing SVD of subject ', subj, '...']), logFile);        
	        [u,s,v] = svd(T, 'econ');
	        s = diag(s);
    	    clear u;
		displayLogItem(sprintf('\t%s', ['Completed computing SVD of subject ', subj]), logFile);
        
        if pctEnergy < 100
			displayLogItem(sprintf('\t%s', ['Projecting subject ', subj, '...']), logFile);        	
			p = findPercentile(s, pctEnergy);
			displayLogItem(sprintf('\t%s', ['Keeping ', num2str(p), ' components']), logFile);
			s = s(1:p);
			v = v(:, 1:p);
			displayLogItem(sprintf('\t%s', ['Completed projecting subject ', subj]), logFile);        	
        end

        % Save out data to directory
        saveOrigProjectedData(subj, v, s);
    end
    displayLogItem('Completed pre-processing the datasets...', logFile);
    clear s v p T;
%% END Data preprocessing
    

%% Determine the subject ordering
	displayLogItem(['Subject ordering:  ', mat2str(cell2mat(subjects'))], logFile);
%% END Determine the subject ordering


%% Determine the starting point, if necessary
    startPass = 0;
    startSubjNum = 1;
    if determineStartPoint == 1
        % If the directory already exists, then try to determine the starting
        % point
        if checkWarpCompleted == 1
            displayLogItem('This atlas construction has completed already.');
            return;
        end
        
        displayLogItem('This atlas has already been partially constructed.  The algorithm will determine the last point and continue from there', logFile);
        foundStart = 0;
        for passNum = 0:numPasses
            for subjNum = 1:length(subjects)
                % Determine if the warp for this subject on this pass exists
                if checkWarpFileExists(subjects{subjNum}, hems, passNum) == 1
                    % Then it exists and is a file
                    continue;
                else
                    % Then file doesn't exist
                    % Break from the loop
                    foundStart = 1;
                    break;
                end
            end

            if foundStart == 1
                break;
            end
        end

        if foundStart == 1
            % Then set the starting point
            startPass = passNum;
            startSubjNum = subjNum;
            displayLogItem(['The algorithm will start on pass #', num2str(startPass), ' and subject #', num2str(startSubjNum), ' (', subjects{startSubjNum}, ') out of ', num2str(length(subjects))], logFile);
        end
    end
    clear determineStartPoint;
%% END Determine the starting point


%% 0th pass
    if startPass == 0
        displayLogItem('Running through pass #0...', logFile);

		% Build up the atlas
		displayLogItem('Building leave-one-out atlas', logFile);
		touched=0;
		nsubjInAtlas = 0;
		for subjNum = 1:nsubj
			if subjNum == startSubjNum
				continue;
			end

			subj = subjects{subjNum};			
			
			if subjNum > startSubjNum
				% Then we must load the original data
				subDataTS = loadOrigProjectedData(subj);
				subDataTS = normalizeDataset(transpose(subDataTS));
				[u,s,v] = svd(subDataTS, 'econ');	
				subData.V = v;
				subData.sigma = diag(s);
				clear u s v subDataTS;

	            saveWarpedProjectedData(subj, subData.V, subData.sigma);
			else
				% Then we must load the warped data
				[subData.V, subData.sigma] = loadWarpedProjectedData(subj);
			end
			
			if touched == 1
				displayLogItem(['Updating the atlas by adding in subject #', num2str(subjNum)], logFile);
				[atlas.V, atlas.sigma] = updateGroupTemplate(atlas.V, atlas.sigma, nsubjInAtlas, subData.V, subData.sigma);
				displayLogItem(['Completed updating the atlas by adding in subject #', num2str(subjNum)], logFile);
			else
				displayLogItem(['Temporarily setting the atlas as subject #', num2str(subjNum)], logFile);
				atlas.V = subData.V;
				atlas.sigma = subData.sigma;
			end
			nsubjInAtlas = nsubjInAtlas + 1;
			touched = 1;
		end
		clear subData;
		displayLogItem('Completed building leave-one-out atlas', logFile);
			
        % Now run through the 0th pass
        for subjNum = startSubjNum:nsubj
            subj = subjects{subjNum};

            % Load the data in
            [subOrig.V, subOrig.sigma] = loadOrigProjectedData(subj);   
            % Set the output file prefix
            alignLogFile = [logDir, 'log_alignment_pass_0_subjNum_', num2str(subjNum), '_subj_', subj, '.txt'];
            % Now run the algorithm
            displayLogItem(['Performing alignment of subject #', num2str(subjNum), ' (', subj, ') on pass #0'], logFile);
            displayLogItem(['A log of this alignment will be saved to ', alignLogFile], logFile);
            
            if strcmp(atlasProjMethod, 'keepAll')
            	warp = connnorm_register(subOrig, atlas, coords, zeros(3, numNodes, opsCoordsType), regularization, maxResolution, alignLogFile);            
            elseif strcmp(atlasProjMethod, 'keepSameAsSubject')
	            ncomponents = min(size(subOrig.V, 2), size(atlas.V, 2));
	            displayLogItem(['Keeping ', num2str(ncomponents), ' of the atlas'], logFile);
    	        atlasReg.V = atlas.V(:, 1:ncomponents);
        	    atlasReg.sigma = atlas.sigma(1:ncomponents);
            	warp = connnorm_register(subOrig, atlasReg, coords, zeros(3, numNodes, opsCoordsType), regularization, maxResolution, alignLogFile);
            	clear atlasReg;
            elseif strcmp(atlasProjMethod, 'keepSpecified')
            	ncomponents = min(size(atlas.V, 2), absMaxAtlasComponents);
            	displayLogItem(['Keeping ', num2str(ncomponents), ' of the atlas'], logFile);
    	        atlasReg.V = atlas.V(:, 1:ncomponents);
        	    atlasReg.sigma = atlas.sigma(1:ncomponents);
            	warp = connnorm_register(subOrig, atlasReg, coords, zeros(3, numNodes, opsCoordsType), regularization, maxResolution, alignLogFile);            
            	clear atlasReg;
            else
				displayLogItem(['ERROR:  Unknown atlasProjMethod:  ', atlasProjMethod], logFile);
				error(['ERROR:  Unknown atlasProjMethod:  ', atlasProjMethod]);
			end            
            displayLogItem(['Completed alignment of subject #', num2str(subjNum), ' (', subj, ') on pass #0'], logFile);

            % Now we need to output the warp file
            saveWarp(warp, subj, 0, 0);

            % Now interpolate the data using the warp
            subWarpedTS = computeInterpOnSphere(diagMultPre(subOrig.sigma, transpose(subOrig.V)), coords.cartCoords, coords.neighbors, warp);
            clear subOrig warp;
            subWarpedTS = normalizeDataset(subWarpedTS);
            % Write out the temp dataset to file
            [u,s,v] = svd(subWarpedTS, 'econ');
            subWarped.V = v;
            subWarped.sigma = diag(s);
            clear u s v subWarpedTS;
            saveWarpedProjectedData(subj, subWarped.V, subWarped.sigma);
            
			% Now update the atlas
            if subjNum < nsubj
            	subRemoveInd = subjNum+1;
				displayLogItem(['Updating group template by adding subject #', num2str(subjNum), ' and removing subject #', num2str(subRemoveInd)], logFile);
				[subRemove.V, subRemove.sigma] = loadWarpedProjectedData(subjects{subRemoveInd});
				[atlas.V, atlas.sigma] = updateGroupTemplate(atlas.V, atlas.sigma, (nsubj-1), subWarped.V, subWarped.sigma, subRemove.V, subRemove.sigma);
				displayLogItem(['Completed updating group template!'], logFile);
				clear subRemove subWarped;			
			end
        end
        
		displayLogItem('Zero correcting the warps, and re-adjusting the atlas...', logFile);
		% Then zero correct the warps
		for subjNum = 1:nsubj
			subj = subjects{subjNum};
			warps{subjNum} = loadWarp(subj, 0, 0);
		end

		warps_zero = computeZeroCorrection(coords.cartCoords, warps);
        for subjNum = 1:nsubj
            subj = subjects{subjNum};
            warp = warps_zero{subjNum};

            saveWarp(warp, subj, 0, 1);
                
            % Re-warp the subject
            [subOrig.V, subOrig.sigma] = loadOrigProjectedData(subj);
            subWarpedTS = computeInterpOnSphere(diagMultPre(subOrig.sigma, transpose(subOrig.V)), coords.cartCoords, coords.neighbors, warp);
            clear warp subOrig;
            subWarpedTS = normalizeDataset(subWarpedTS);
            [u,s,v] = svd(subWarpedTS, 'econ');
            subWarped.sigma = diag(s);
            subWarped.V = v;
            clear u s v subWarpedTS;
            % Write out the warped dataset
            saveWarpedProjectedData(subj, subWarped.V, subWarped.sigma);
             
            % Now we also need to adjust the atlas...
            if subjNum == 1
            	% Do nothing -- we don't want the first subject in the atlas right now since we will be starting a new pass
            elseif subjNum == 2
            	displayLogItem('Initially setting subject #2 as the atlas', logFile);
            	atlas.V = subWarped.V;
            	atlas.sigma = subWarped.sigma;
            else
            	displayLogItem(['Updating the group template by adding in subject #', num2str(subjNum)], logFile);
            	[atlas.V, atlas.sigma] = updateGroupTemplate(atlas.V, atlas.sigma, (subjNum - 2), subWarped.V, subWarped.sigma);
            	displayLogItem(['Completed updating the group template by adding in subject #', num2str(subjNum)], logFile);
            end
            clear subWarped;
            displayLogItem('Completed zero correcting the warps, and re-adjusting the atlas', logFile);
        end        
        
        % Restart start values for next set of passes
        startPass = 1;
        startSubjNum = 1;
        displayLogItem('Completed running through pass #0', logFile);
    else
        % Then we need to develop the atlas
        % Also need to load in the warps themselves
        % Run through all subjects
        displayLogItem('Re-creating atlas from previous run...', logFile);
        [atlas.V, atlas.sigma] = loadWarpedProjectedData(subjects{2});
        for subjNum = 3:nsubj
            [VA, sa] = loadWarpedProjectedData(subjects{subjNum});
            [atlas.V, atlas.sigma] = updateGroupTemplate(atlas.V, atlas.sigma, (subjNum-2), VA, sa);
        end
        clear VA sa;
        displayLogItem('Completed re-creating atlas from previous run', logFile);
    end
%% END 0th pass


%% Remaining passes
    % Now run through the remaining passes
    for passNum = startPass:numPasses
        displayLogItem(['Running through pass #', num2str(passNum), '...'], logFile);
        for subjNum = startSubjNum:nsubj
            subj = subjects{subjNum};

			% Load in the original projected data for this subject
			[subOrig.V, subOrig.sigma] = loadOrigProjectedData(subj);

            % Set the output file prefix
            alignLogFile = [logDir, 'log_alignment_pass_', num2str(passNum), '_subjNum_', num2str(subjNum), '_subj_', subj, '.txt'];
            % Load the current warp file 
            if passNum > 1
                warp = loadWarp(subj, passNum-1, 1);
            else
                warp = loadWarp(subj, passNum-1, 0);
            end
            % Now run the algorithm            
            displayLogItem(['Performing alignment of subject #', num2str(subjNum), ' (', subj, ') on pass #', num2str(passNum)], logFile);
            displayLogItem(['A log of this alignment will be saved to ', alignLogFile], logFile);

            if strcmp(atlasProjMethod, 'keepAll')
	            warp = connnorm_register(subOrig, atlas, coords, warp, regularization, maxResolution, alignLogFile);
            elseif strcmp(atlasProjMethod, 'keepSameAsSubject')
	            ncomponents = min(size(subOrig.V, 2), size(atlas.V, 2));
	            displayLogItem(['Keeping ', num2str(ncomponents), ' of the atlas'], logFile);
    	        atlasReg.V = atlas.V(:, 1:ncomponents);
        	    atlasReg.sigma = atlas.sigma(1:ncomponents);
	            warp = connnorm_register(subOrig, atlasReg, coords, warp, regularization, maxResolution, alignLogFile);
	            clear atlasReg;
	        elseif strcmp(atlasProjMethod, 'keepSpecified')
            	ncomponents = min(size(atlas.V, 2), absMaxAtlasComponents);
	            displayLogItem(['Keeping ', num2str(ncomponents), ' of the atlas'], logFile);
    	        atlasReg.V = atlas.V(:, 1:ncomponents);
        	    atlasReg.sigma = atlas.sigma(1:ncomponents);
	            warp = connnorm_register(subOrig, atlasReg, coords, warp, regularization, maxResolution, alignLogFile);	        
	            clear atlasReg;
            else
				displayLogItem(['ERROR:  Unknown atlasProjMethod:  ', atlasProjMethod], logFile);
				error(['ERROR:  Unknown atlasProjMethod:  ', atlasProjMethod]);
			end     

            % Now we need to output the warp file
            saveWarp(warp, subj, passNum, 0);
            % Now interpolate the data using the warp
            subWarpedTS = computeInterpOnSphere(diagMultPre(subOrig.sigma, transpose(subOrig.V)), coords.cartCoords, coords.neighbors, warp);
            clear warp subOrig;
            subWarpedTS = normalizeDataset(subWarpedTS);
            [u,s,v] = svd(subWarpedTS, 'econ');
            subWarped.sigma = diag(s);
            subWarped.V = v;
            clear u s v subWarpedTS;

            % Write out the temp dataset to file
            saveWarpedProjectedData(subj, subWarped.V, subWarped.sigma);

            % Now update the atlas
            if subjNum < nsubj
            	[subNext.V, subNext.sigma] = loadWarpedProjectedData(subjects{subjNum+1});
            	[atlas.V, atlas.sigma] = updateGroupTemplate(atlas.V, atlas.sigma, (nsubj-1), subWarped.V, subWarped.sigma, subNext.V, subNext.sigma);
            	clear subNext;
            end
        end

        if passNum > 0
            displayLogItem('Zero correcting the warps, and re-adjusting the atlas...', logFile);
            % Then zero correct the warps
            for subjNum = 1:nsubj
                subj = subjects{subjNum};
                warps{subjNum} = loadWarp(subj, passNum, 0);
            end

            warps_zero = computeZeroCorrection(coords.cartCoords, warps);
            for subjNum = 1:nsubj
                subj = subjects{subjNum};
                warp = warps_zero{subjNum};

                saveWarp(warp, subj, passNum, 1);
                
                % Re-warp the subject
                [subOrig.V, subOrig.sigma] = loadOrigProjectedData(subj);
                subWarpedTS = computeInterpOnSphere(diagMultPre(subOrig.sigma, transpose(subOrig.V)), coords.cartCoords, coords.neighbors, warp);
                clear warp subOrig;
                subWarpedTS = normalizeDataset(subWarpedTS);
                [u,s,v] = svd(subWarpedTS, 'econ');
                subWarped.sigma = diag(s);
                subWarped.V = v;
                clear u s v subWarpedTS;
                % Write out the warped dataset
                saveWarpedProjectedData(subj, subWarped.V, subWarped.sigma);
                
                if passNum < numPasses
                	% Then we also need to adjust the atlas
                	if subjNum == 1
                		% Do nothing
                	elseif subjNum == 2
                		atlas.V = subWarped.V;
                		atlas.sigma = subWarped.sigma;
                	else
                		[atlas.V, atlas.sigma] = updateGroupTemplate(atlas.V, atlas.sigma, (subjNum - 2), subWarped.V, subWarped.sigma);
                	end
                end                
            end        
            clear subWarped;
            displayLogItem('Completed zero correcting the warps, and re-adjusting the atlas', logFile);
        end
        
        % Reset the subject starting number
        startSubjNum = 1;
        displayLogItem(['Completed running through pass #', num2str(passNum)], logFile);
    end
%% END Remaining passes

    % Save final warps
    for subjNum = 1:nsubj
        subj = subjects{subjNum};
        saveFinalWarp(subj);
    end
    
    displayLogItem('Completed alignment', logFile);
    








%% FUNCTION getOrigTimeSeriesFilename
    function fName = getOrigTimeSeriesFilename(subj, hem)
        fName = [outDir, subj, '_', hem, '_', num2str(mmRes), 'mm_fwhm', num2str(fwhm), '_', experiment, '_on_', fs_surf, '.niml.dset'];
    end
%% END FUNCTION getOrigTimeSeriesFilename


%% FUNCTION getOrigProjectedDataFilenames
    function [V_fName, s_fName] = getOrigProjectedDataFilenames(subj)
        V_fName = [tmpDir, 'V_', subj, '_orig.bin'];
        s_fName = [tmpDir, 's_', subj, '_orig.bin'];
    end
%% END FUNCTION getOrigProjectedDataFilenames


%% FUNCTION getWarpedProjectedDataFilenames
    function [V_fName, s_fName] = getWarpedProjectedDataFilenames(subj)
        V_fName = [tmpDir, 'V_', subj, '_warped.bin'];
        s_fName = [tmpDir, 's_', subj, '_warped.bin'];
    end
%% END FUNCTION getTmpTimeseriesFilenames


%% FUNCTION loadOrigTimeseries
    function T = loadOrigTimeseries(subj)
%        T = zeros(tLen, nhem*numNodesPerHem, opsDataType);

        n = 1;
        for hemNum = 1:nhem
            hem = hems{hemNum};
            T(:, n:(n+numNodesPerHem-1)) = read_niml_binary(getOrigTimeSeriesFilename(subj, hem), opsDataType, saveDataType);

            n = n + numNodesPerHem;
        end
        T = cast(T(1:tLen, :), opsDataType);
    end
%% END FUNCTION loadOrigTimeseries


%% FUNCTION loadOrigProjectedData
    function [V, s] = loadOrigProjectedData(subj)
        [V_fName, s_fName] = getOrigProjectedDataFilenames(subj);
        [V, s] = loadProjectedData(subj, V_fName, s_fName);
        
        if nargout == 1
            V = diagMultPost(V, s);
            s = [];
        end
    end
%% END FUNCTION loadOrigProjectedData


%% FUNCTION loadWarpedProjectedData
    function [V, s] = loadWarpedProjectedData(subj)
        [V_fName, s_fName] = getWarpedProjectedDataFilenames(subj);
        [V, s] = loadProjectedData(subj, V_fName, s_fName);
    end
%% END FUNCTION loadWarpedProjectedData


%% FUNCTION loadProjectedData
    function [V, s] = loadProjectedData(subj, V_fName, s_fName)
        d = dir(s_fName);
        numEigs = (d.bytes)/saveDataBytes;
        
        fp = fopen(s_fName, 'rb');
        s = fread(fp, numEigs, [saveDataType, '=>', opsDataType]);
        fclose(fp);
        
        fp = fopen(V_fName, 'rb');
        V = fread(fp, [numNodes, numEigs], [saveDataType, '=>', opsDataType]);
        fclose(fp);        
    end
%% END FUNCTION loadProjectedData


%% FUNCTION saveOrigProjectedData
    function saveOrigProjectedData(subj, V, s)
        [V_fName, s_fName] = getOrigProjectedDataFilenames(subj);
        saveProjectedData(V, s, V_fName, s_fName);
    end
%% END FUNCTION saveOrigProjectedData


%% FUNCTION saveWarpedProjectedData
    function saveWarpedProjectedData(subj, V, s)
        [V_fName, s_fName] = getWarpedProjectedDataFilenames(subj);
        saveProjectedData(V, s, V_fName, s_fName);
    end
%% END FUNCTION saveWarpedProjectedData


%% FUNCTION saveProjectedData
    function saveProjectedData(V, s, V_fName, s_fName)
        fp = fopen(V_fName, 'wb');
        fwrite(fp, V, saveDataType);
        fclose(fp);
        
        fp = fopen(s_fName, 'wb');
        fwrite(fp, s, saveDataType);
        fclose(fp);        
    end
%% END FUNCTION saveProjectedData


%% FUNCTION checkWarpFileExists
    function e = checkWarpFileExists(subj, hems, passNum)
        e = 1;
        for hemNum = 1:nhem
            hem = hems{hemNum};

            ex = exist(getWarpFilename(subj, hem, passNum, 0));
            if ex == 0
                e = 0;
                break;
            elseif ex == 2
                % Then file exists
            else
                error('Encountered strange exist result');
            end
        end
    end
%% END FUNCTION checkWarpFilesExists


%% FUNCTION getWarpFilename
    function warpFilename = getWarpFilename(subj, hem, passNum, zeroCorrected)
        if nargin == 3
            zeroCorrected = 0;
        end

        corrStr = '';
        if zeroCorrected == 1
            corrStr = '_corrected';
        end

        warpFilename = [tmpDir, 'warp_', subj, '_', hem, '_', num2str(passNum), corrStr, '.bin'];
    end
%% END FUNCTION getWarpFilename


%% FUNCTION loadWarp
    function warp = loadWarp(subj, passNum, zeroCorrected)
        warp = zeros(3, nhem*numNodesPerHem, opsCoordsType);

        n = 1;
        for hemNum = 1:nhem
            hem = hems{hemNum};

            fp = fopen(getWarpFilename(subj, hem, passNum, zeroCorrected), 'rb');
            warp(:, n:(n+numNodesPerHem-1)) = fread(fp, [3, numNodesPerHem], [saveCoordsType, '=>', opsCoordsType]);
            fclose(fp);

            n = n + numNodesPerHem;
        end
    end
%% END FUNCTION loadWarp


%% FUNCTION saveWarp
    function saveWarp(warp, subj, passNum, zeroCorrected)
        n = 1;
        for hemNum = 1:nhem
            hem = hems{hemNum};

            fp = fopen(getWarpFilename(subj, hem, passNum, zeroCorrected), 'wb');
            fwrite(fp, warp(:, n:(n+numNodesPerHem-1)), saveCoordsType);
            fclose(fp);

            n = n + numNodesPerHem;
        end
    end
%% END FUNCTION saveWarp


%% FUNCTION saveFinalWarp
    function saveFinalWarp(subj)
        warp = loadWarp(subj, numPasses, 1);
        
        n = 1;
        for hemNum = 1:nhem
            hem = hems{hemNum};
            hemDir = hemOutDirs{hemNum};
            
            fp = fopen([hemDir, 'warp_', subj, '_', hem, '_final.bin'], 'wb');
            wh = warp(:, n:(n+numNodesPerHem-1));
            fwrite(fp, wh, saveCoordsType);
            fclose(fp);
            
            fp = fopen([hemDir, 'standard', num2str(mmRes), 'mm_', subj, '_', hem, '.connnorm.asc'], 'w');
            fprintf(fp, '#\n');            
            fprintf(fp, '%d %d\n', numNodesPerHem, numTrisPerHem);
            fprintf(fp, '%f %f %f 0\n', 100*(coords.cartCoords(:, n:(n+numNodesPerHem-1)) + wh));
            fprintf(fp, '%d %d %d 0\n', coords.triangles(:, 1:numTrisPerHem));
            fclose(fp);
            
            n = n + numNodesPerHem;
        end
    end
%% END FUNCTION saveFinalWarp


%% FUNCTION checkWarpCompleted
    function completed = checkWarpCompleted
        completed = 1;
        for subjNum = 1:nsubj
            for hemNum = 1:nhem
                hem = hems{hemNum};
                hemDir = hemOutDirs{hemNum};
                
                ex = exist([hemDir, 'warp_', subj, '_', hem, '_final.bin']);
                if ex == 0
                    completed = 0;
                    return;
                elseif ex == 2
                    % Then file exists
                else
                    error('Strange exist number');
                end
            end
        end
    end
%% END FUNCTION checkWarpCompleted

end
%% END connnorm