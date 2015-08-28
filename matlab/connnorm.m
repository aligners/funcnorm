function warps = connnorm(subjects, hems, experiment, lambda_areal, lambda_metric, varargin)
% FUNCTION warps = connnorm(subjects, hems, experiment, lambda_areal, lambda_metric, varargin)
% *** INPUT ARGUMENTS ***
%   subjects:  a cell array of strings specifying the subject names
%       Ex:  {'cb', 'dm', 'hj', 'kd', 'kl', 'mh', 'ph', 'rb', 'se', 'sm'}
%   hems:  a cell array of strings specifying the hemispheres to align
%       Ex:  {'lh', 'rh'}
%	experiment:  name of the functional experiment (e.g., 'movie5min_P1')
%
%	connnorm(subjects, hems, experiment, 'NAME', VALUE, ...) allows you to
%	specify additional options
%		in name/value pairs as follows:
%
%		'outDir' - the specified output directory (DEFAULT:  '.')
%		'alignDir' - alignment directory name - will be relative to outDir (DEFAULT:  'alignment_%s')
%		'alignment_continue' - a flag to continue a partially completed alignment (DEFAULT:  0)
%			If alignment_continue is set to 0 and the specified align directory exists already, data will be overwritten
%			If alignment_continue is set to 1 and the specified align directory exists already, then will try to continue where it left off
%		'logFile' - filename to store log items (relative to sprintf('%s/log', alignDir)) (DEFAULT:  'connnorm_log.txt')
%		'mmRes' - the resolution of the cortical mesh (DEFAULT:  2)
%		'fwhm' - the fwhm for SurfSmooth (DEFAULT:  0)
%		'fs_surf' - the freesurfer surface file (e.g., 'sphere', or 'sphere.reg') (DEFAULT:  'sphere.reg')
%		'pctEnergy' - the percentage of fMRI time-series energy to preserve in the PCA dimensionality
%			reduction step (DEFAULT:  100)
%		'maxResolution' - multiplier specifying the maximum resolution for multi-resolution alignment.
%			The alignment begins at a resolution of maxResolution*mmRes, where mmRes is specified as the
%			initial cortical mesh resolution (DEFAULT:  3)
%		'numPasses' - the number of alignment passes through the subjects (DEFAULT:  3)
%		'connnormDir' - string specifying the full path of connnorm_release (DEFAULT:  '.')
%		'startFromZScoredData' - specify whether or not to z-score the original data (DEFAULT:  1)
%		'runFuncnormPrepare' - flag specifying whether or not to run funcnorm_prepare script (DEFAULT:  0)
%		'tWin' - a vector of time indices to include in the connectivity alignment algorithm
%		'useMultithreading' - a flag to specify whether or not to enable MATLAB's multithreading option (DEFAULT:  1)
%		'useGroupSVDRows' - a flag specifying how to compress the fMRI datasets.  If this flag is set to true, then a common set of singular vectors (that define the subspace) are computed with respect to all the subjects, and then each subject's data is projected into this subspace
%			Otherwise, if useGroupSVD is set to 0, then a SVD is run on each individual subject.
%			useGroupSVDRows=1 is recommended if pctEnergy is low, because it minimizes the lack of overlap of subspaces that may occur if each SVD is performed individually.
%       'useGroupSVDCols' - a flag specifying how to compress the fMRI datasets.
%       'initWarp' - a cell array specifying an initial set of warps to use
%           initWarp{j} would correspond to the warp for subjects{j}  
%
%	If runFuncnormPrepare is set, then the following options may also be specified:
%
%		'volDataDir' - the directory location to the volume data
%			Ex:  '/data/subjects_in/sample_data/movie' (REQUIRED)
%		'inpSubjDir' - input subjects directory 
%			Ex:  '/data/subjects_in/sample_data/SUBJECTS_DIR' (REQUIRED)
%		'saveSubjDir' - directory to save generated surfaces and intermediate results (REQUIRED)
%		'suffix' - suffix to volume names (DEFAULT:  '+orig')
%		'extensions' - the extension to the anatomical and functional volumes (DEFAULT:  'BRIK.gz')
%		'preprocess_movie' - a flag (0 or 1) specifying whether to perform movie preprocessing 
%			on the functional data (DEFAULT:  0)
%
% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

	% Clear all breakpoints in case there was some earlier debugging going on...
	dbclear all;

%% Initialize all variables to defaults
	outDir = '.';
	alignDir = sprintf('alignment_%s', experiment);
	logFile = 'connnorm_log.txt';
	startFromZScoredData = 1;
	connnormDir = '.';
	mmRes = 2;
	fwhm = 0;
	fs_surf = 'sphere.reg';
	pctEnergy = 100;
	maxResolution = 3;
	numPasses = 3;
	firstPassAtlasMethod = 'average'; % Either 'build_up' or 'average'
	runFuncnormPrepare = 0;
	alignment_continue = 0;
	suffix = '+orig';
	extensions = 'BRIK.gz';
	preprocess_movie = 0;
    saveDataBytes = 4; % Number of bytes of datatype to store data (e.g. float32 = 4 bytes)
    saveDataType = 'single';
    opsDataBytes = 8; % Number of bytes to perform all operations in (e.g. single precision = 4 bytes)
    opsDataType = 'single';
    opsCoordsType = 'double'; % Format for coords structures
    saveCoordsType = 'single';
    useTWin = 0;
    useGroupSVDRows = 0;
    useGroupSVDCols = 0;
    useMultithreading = 1;
    useInitWarp = 0;
%% Parse property name/value pair input arguments
	i = 1;
	while i <= length(varargin)
		propName = varargin{i};
		switch propName
			case 'connnormDir'
				connnormDir = varargin{i+1};
				i = i+2;
			case 'startFromZScoredData'
				startFromZScoredData = varargin{i+1};
				i = i+2;
			case 'tWin'
				tWin = varargin{i+1};
				useTWin = 1;
				i = i+2;
				if length(tWin) == 1
					if tWin == -1
						useTWin = 0;
					else
						useTWin = 1;
					end
				end
			case 'useGroupSVDRows'
				useGroupSVDRows = varargin{i+1};
				i = i+2;
            case 'useGroupSVDCols'
                useGroupSVDCols = varargin{i+1};
                i = i+2;
			case 'alignment_continue'
				alignment_continue = varargin{i+1};
				i = i+2;
			case 'outDir'
				outDir = varargin{i+1};
				i = i+2;
			case 'alignDir'
				alignDir = varargin{i+1};
				i = i+2;
			case 'logFile'
				logFile = varargin{i+1};
				i = i+2;
			case 'mmRes'
				mmRes = varargin{i+1};
				i = i+2;
			case 'fwhm'
				fwhm = varargin{i+1};
				i = i+2;
			case 'fs_surf'
				fs_surf = varargin{i+1};
				i = i+2;
			case 'pctEnergy'
				pctEnergy = varargin{i+1};
				i = i+2;
			case 'maxResolution'
				maxResolution = varargin{i+1};
				i = i+2;
			case 'numPasses'
				numPasses = varargin{i+1};
				i = i+2;
			case 'firstPassAtlasMethod'
				firstPassAtlasMethod = varargin{i+1};
				i = i+2;
			case 'runFuncnormPrepare'
				runFuncnormPrepare = varargin{i+1};
				i = i+2;
			case 'volDataDir'
				volDataDir = varargin{i+1};
				i = i+2;
			case 'suffix'
				suffix = varargin{i+1};
				i = i+2;
			case 'extensions'
				extensions = varargin{i+1};
				i = i+2;
			case 'inpSubjDir'
				inpSubjDir = varargin{i+1};
				i = i+2;
			case 'saveSubjDir'
				saveSubjDir = varargin{i+1};
				i = i+2;
			case 'preprocess_movie'
				preprocess_movie = varargin{i+1};
				i = i+2;
			case 'useMultithreading'
				useMultithreading = varargin{i+1};
				i = i+2;
            case 'initWarp'
                initWarp = varargin{i+1};
                useInitWarp = 1;
                i = i+2;
			otherwise
				error(sprintf('Unknown property name [%s]', propName));
		end	
    end	
    if useGroupSVDRows == 1 && useGroupSVDCols == 1
        error('Cannot set both useGroupSVDRows and useGroupSVDCols flags');
    end
%% END Parse property name/value pair input arguments
        
    if isdir([connnormDir, '/mex_files']) == 0
    	error(['Cannot find mex_files directory (looking here:  ', connnormDir, '). Consider changing the connnormDir setting.']);
    end
    if isdir([connnormDir, '/m_files']) == 0
    	error(['Cannot find m_files directory (looking here:  ', connnormDir, '). Consider changing the connnormDir setting.']);
    end
    addpath([connnormDir, '/mex_files/']);
    addpath([connnormDir, '/m_files/']);

%% SETUP output directory structure
	if isdir(outDir) == 0; mkdir(outDir); end;
	alignDirOrig = alignDir;
	alignDir = sprintf('%s/%s', outDir, alignDir);
	continueFromWhereLeftOff = 0;
	if isdir(alignDir) == 0; 
		mkdir(alignDir); 
	else
		if alignment_continue == 0
			warning(['Alignment directory already exists.  Contents will be overwritten (in the future, to continue where you left off, set the alignment_continue flag)']);
			system(sprintf('rm -rf %s/*', alignDir));
		else
			continueFromWhereLeftOff = 1;
		end
	end;
	% confDir is where all constants are saved
	confDir = sprintf('%s/conf', alignDir); if isdir(confDir) == 0; mkdir(confDir); end;
	% logDir is where all logs are stored
	logDir = sprintf('%s/log', alignDir);if isdir(logDir) == 0; mkdir(logDir); end;
	% tmpDir is a temporary directory for this alignment
	tmpDir = sprintf('%s/tmp', alignDir);if isdir(tmpDir) == 0; mkdir(tmpDir); end;
	% tmpDir2 is a second temporary directory for this alignment
	tmpDir2 = sprintf('%s/tmp2/%s', alignDir, alignDirOrig); if isdir(tmpDir2) == 0; mkdir(tmpDir2); end;
	% warpsDir is where all final warps are stored
	warpsDir = sprintf('%s/warps', alignDir);if isdir(warpsDir) == 0; mkdir(warpsDir); end;
	% logFile is the name of the main algorithm log file
	logFile = sprintf('%s/%s', logDir, logFile);
%% END Setup output directory structure
    
    if continueFromWhereLeftOff == 1
    	% Then we need to check if any configuration settings have changed in the .conf file
    	% If so, then we must error
		fp = fopen(sprintf('%s/alignment_tmp.conf', confDir), 'w');
		fprintf(fp, 'subjects:%s\n', mat2str(cell2mat(subjects')));
		fprintf(fp, 'hems:%s\n', mat2str(cell2mat(hems')));
		fprintf(fp, 'experiment:%s\n', experiment);
		fprintf(fp, 'lambda_areal:%f\n', lambda_areal);
		fprintf(fp, 'lambda_metric:%f\n', lambda_metric);
		fprintf(fp, 'logFile:%s\n', logFile);
		fprintf(fp, 'mmRes:%d\n', mmRes);
		fprintf(fp, 'fwhm:%f\n', fwhm);
		fprintf(fp, 'fs_surf:%s\n', fs_surf);
		fprintf(fp, 'pctEnergy:%f\n', pctEnergy);
		fprintf(fp, 'maxResolution:%d\n', maxResolution);
		fprintf(fp, 'numPasses:%d\n', numPasses);
		fprintf(fp, 'firstPassAtlasMethod:%s\n', firstPassAtlasMethod);
		fprintf(fp, 'runFuncnormPrepare:%d\n', runFuncnormPrepare);
		fprintf(fp, 'suffix:%s\n', suffix);
		fprintf(fp, 'extensions:%s\n', extensions);
		fprintf(fp, 'preprocess_movie:%d\n', preprocess_movie);
		fprintf(fp, 'saveDataType:%s\n', saveDataType);
		fprintf(fp, 'opsDataType:%s\n', opsDataType);
		fprintf(fp, 'saveCoordsType:%s\n', saveCoordsType);
		fprintf(fp, 'opsCoordsType:%s\n', opsCoordsType);
		fprintf(fp, 'startFromZScoredData:%d\n', startFromZScoredData);
		if useTWin == 1
			fprintf(fp, 'time window:%s\n', mat2str(tWin));
		end
		if useGroupSVDRows == 1
			fprintf(fp, 'Using group SVD Rows\n');
		end
		fclose(fp);	
		
		[status,result] = system(sprintf('diff %s/alignment.conf %s/alignment_tmp.conf', confDir, confDir));
		if status >= 2
			error(sprintf('Couldn''t find original conf file in %s directory. Cannot continue from previous alignment', confDir));
		elseif status == 1
			% Then there are differences
			error(sprintf('Parameter settings have changed from the previous alignment. Please check the configuration file in %s directory.', confDir));
		else
			% Everything is OK...
			displayLogItem('Continuing connnorm algorithm from where it left off...', logFile);
		end
			
		system(sprintf('rm -f %s/alignment_tmp.conf', confDir));
    end
    
    if startFromZScoredData == 0
        error('startFromZScoredData == 0 is no longer supported');
    end

%% WRITE out configuration settings to .conf file
	fp = fopen(sprintf('%s/alignment.conf', confDir), 'w');
	fprintf(fp, 'subjects:%s\n', mat2str(cell2mat(subjects')));
	fprintf(fp, 'hems:%s\n', mat2str(cell2mat(hems')));
	fprintf(fp, 'experiment:%s\n', experiment);
	fprintf(fp, 'lambda_areal:%f\n', lambda_areal);
	fprintf(fp, 'lambda_metric:%f\n', lambda_metric);
	fprintf(fp, 'logFile:%s\n', logFile);
	fprintf(fp, 'mmRes:%d\n', mmRes);
	fprintf(fp, 'fwhm:%f\n', fwhm);
	fprintf(fp, 'fs_surf:%s\n', fs_surf);
	fprintf(fp, 'pctEnergy:%f\n', pctEnergy);
	fprintf(fp, 'maxResolution:%d\n', maxResolution);
	fprintf(fp, 'numPasses:%d\n', numPasses);
	fprintf(fp, 'firstPassAtlasMethod:%s\n', firstPassAtlasMethod);
	fprintf(fp, 'runFuncnormPrepare:%d\n', runFuncnormPrepare);
	fprintf(fp, 'suffix:%s\n', suffix);
	fprintf(fp, 'extensions:%s\n', extensions);
	fprintf(fp, 'preprocess_movie:%d\n', preprocess_movie);
	fprintf(fp, 'saveDataType:%s\n', saveDataType);
	fprintf(fp, 'opsDataType:%s\n', opsDataType);
	fprintf(fp, 'saveCoordsType:%s\n', saveCoordsType);
	fprintf(fp, 'opsCoordsType:%s\n', opsCoordsType);
	fprintf(fp, 'startFromZScoredData:%d\n', startFromZScoredData);
	if useTWin == 1
		fprintf(fp, 'time window:%s\n', mat2str(tWin));
	end
	if useGroupSVDRows == 1
		fprintf(fp, 'Using group SVD Rows\n');
    end
    if useGroupSVDCols == 1
        fprintf(fp, 'Using group SVD Cols\n');
    end
	fclose(fp);
%% END WRITE out configuration settings to .conf file


    nsubj = length(subjects); % Number of subjects
    nhem = length(hems); % Number of hemispheres
    for hemNum = 1:nhem
        hem = hems{hemNum};

        hemOutDirs{hemNum} = sprintf('%s/%s', warpsDir, hem);
        if isdir(hemOutDirs{hemNum}) == 0
            mkdir(hemOutDirs{hemNum});
        end
    end


    displayLogItem('Beginning connnorm algorithm...', logFile); 
    
    % Check if multithreading is on...
    mCT = maxNumCompThreads;
    if mCT == 1 && useMultithreading == 1
    	displayLogItem('Multithreading is currently set to off.  Currently resetting multithreading to on.', logFile);
 		mCT = maxNumCompThreads('automatic');   	
%    	cTRep = input('Multithreading is currently set to off.  It is strongly recommended that it be turned on.  Would you like to turn it on? Y/N [Y]:  ', 's');
%    	if isempty(cTRep)
%    		cTRep = 'Y';
%    	end
    	
%    	if cTRep == 'Y';
%    		mCT = maxNumCompThreads('automatic');
%		end
	elseif mCT > 1 && useMultithreading == 0
		displayLogItem('Multithreading is currently turned on, but you have selected not to use multithreading.  Turning off multithreading now.', logFile);
		mCT = maxNumCompThreads(1);
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
    	displayLogItem(sprintf('\t%s', ['First pass atlas method:  ', firstPassAtlasMethod]), logFile);
    	displayLogItem(sprintf('\t%s', ['startFromZScoredData:  ', num2str(startFromZScoredData)]), logFile);
    	if useTWin == 1
    		displayLogItem(sprintf('\t%s', ['time window:  ', mat2str(tWin)]), logFile);
    	end
		if useGroupSVDRows == 1
			displayLogItem(sprintf('\t%s', ['Using group SVD Rows']), logFile);
        end
        if useGroupSVDCols == 1
            displayLogItem(sprintf('\t%s', ['Using group SVD Cols']), logFile);
        end
        

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
    surfFile = sprintf('%s/standard%dmm_%s.asc', outDir, mmRes, fs_surf);
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
    tLenFull = getNumTimepoints(getOrigTimeSeriesFilename(subjects{1}, hems{1}));
	if useTWin == 1
		tLen = min(length(tWin), tLenFull);
	else
		tLen = tLenFull;
	end
	displayLogItem(['Number of nodes per hemisphere:  ', num2str(numNodesPerHem)], logFile);
	displayLogItem(['Number of triangles per hemisphere:  ', num2str(numTrisPerHem)], logFile);
	displayLogItem(['Number of timepoints:  ', num2str(tLen)], logFile);
	displayLogItem(['Number of hemispheres:  ', num2str(nhem)], logFile);
	
%% Setup regularization
    regularization.mode = 'metric_and_areal_distortion';
    regularization.lambda_areal = numNodes*lambda_areal;
    regularization.lambda_metric = numNodes*lambda_metric;
%% END Setup regularization


    if useInitWarp == 0
        for subjNum = 1:nsubj
            initWarp{subjNum} = zeros(3, numNodes, opsCoordsType);
        end
    end

%% Data preprocessing
    displayLogItem('Pre-processing the datasets...', logFile);    
    foundDatas = zeros(1, nsubj);
    for subjNum = 1:nsubj
        subj = subjects{subjNum};

	    if continueFromWhereLeftOff == 1
	    	displayLogItem(['Checking if original projected data for subject ', subj, ' exists...'], logFile);
    		% Then first check to see if the file exists
    		if checkIfOrigProjectedDataExists(subj)	== 1
    			displayLogItem(['Found original projected data for subject ', subj], logFile);
    			foundDatas(subjNum) = 1;
    		end
    		displayLogItem(['Couldn''t find original projected data for subject ', subj], logFile);
    	end
    end
    
    if (useGroupSVDRows == 1) || (useGroupSVDCols == 1)
    	if sum(foundDatas) < nsubj
    		% First step is to compute the group SVD of all of the subjects
    		displayLogItem(sprintf('\tBuilding the group SVD of all of the subjects...'), logFile);
			for subjNum = 1:nsubj
				subj = subjects{subjNum};
				
   				T = loadOrigTimeseries(subj, startFromZScoredData);
                if useGroupSVDRows == 1
                    % Don't do anything
                elseif useGroupSVDCols == 1
                    % Take the transpose
                    T = T';
                else
                    error('What happened here?');
                end
                    
				% Here's we'll keep 95% of the energy of the datasets when building the group SVD
				displayLogItem(sprintf('\t%s', ['Computing SVD of subject ', subj, '...']), logFile);        
		        	[u,s,v] = svd(T, 'econ');
			        s = diag(s);
				displayLogItem(sprintf('\t%s', ['Completed computing SVD of subject ', subj]), logFile);
        
				displayLogItem(sprintf('\t%s', ['Projecting subject ', subj, '...']), logFile);        	
				p = findPercentile(s, 95);
				displayLogItem(sprintf('\t%s', ['Keeping ', num2str(p), ' components']), logFile);
				s = s(1:p);
				v = v(:, 1:p);
				displayLogItem(sprintf('\t%s', ['Completed projecting subject ', subj]), 	logFile);        	
        
        		if subjNum == 1
        			groupV = v;
        			groupS = s;
        		else
        			[groupV, groupS, hasErrored] = updateGroupTemplate(groupV, groupS, (subjNum-1), v, s);
        			if hasErrored == 1
        				whos
        				error('Update group template didn''t have enough memory to complete');
        			end
        		end
        	end
        	displayLogItem(sprintf('\tCompleted building the group SVD of all of the subjects!'), logFile);
			save([tmpDir2, '/groupS.mat'], 'groupS');
%			save([tmpDir2, '/groupV.mat'], 'groupV');
        	
        	if pctEnergy < 100
        		displayLogItem(sprintf('\t%s', 'Projecting the group SVD...'), logFile);
        		pGroup = findPercentile(groupS, pctEnergy);
        		displayLogItem(sprintf('\t%s', ['Keeping ', num2str(pGroup), ' components']), logFile);
        		groupS = groupS(1:pGroup);
        		groupV = groupV(:, 1:pGroup);
        	else
        		pGroup = length(groupS);
        	end
        	
			% Now we must go through each subject and project its data into the group components
			for subjNum = 1:nsubj
				subj = subjects{subjNum};
				
				displayLogItem(sprintf('\tProjecting data from subject %s onto the group SVD', subj), logFile);
				T = loadOrigTimeseries(subj, startFromZScoredData);
                if useGroupSVDRows == 1
    				T = (T*groupV)*groupV';
                elseif useGroupSVDCols == 1
                    T = groupV'*T;
                else
                    error('What happened here?');
                end
				if startFromZScoredData == 1
					T = normalizeDataset(T);
				end
				[u,s,v] = svd(T,'econ');
				s = diag(s(1:pGroup,1:pGroup));
				v = v(:, 1:pGroup);

	    	    % Save out data to directory
    	    	saveOrigProjectedData(subj, v, s);
			end
			clear groupS groupV u T s v pGroup;
    	end
    else
	    for subjNum = 1:nsubj
    		subj = subjects{subjNum};
    	
			if foundDatas(subjNum) == 1; continue; end;
			
        	% First, load the data and center it
	        T = loadOrigTimeseries(subj, startFromZScoredData);
        
    	    % Compute the singular vectors/values of the time-series matrix
			displayLogItem(sprintf('\t%s', ['Computing SVD of subject ', subj, '...']), logFile);        
	        	[u,s,v] = svd(T, 'econ');
		        s = diag(s);
			displayLogItem(sprintf('\t%s', ['Completed computing SVD of subject ', subj]), logFile);
        
        	if pctEnergy < 100
				displayLogItem(sprintf('\t%s', ['Projecting subject ', subj, '...']), logFile);        	
				p = findPercentile(s, pctEnergy);
				displayLogItem(sprintf('\t%s', ['Keeping ', num2str(p), ' components']), logFile);
				s = s(1:p);
				v = v(:, 1:p);
				displayLogItem(sprintf('\t%s', ['Completed projecting subject ', subj]), logFile);        	
    	    end
        
        	if startFromZScoredData == 1
        		displayLogItem(sprintf('\tRe-normalizing dataset of subject %s...', subj), logFile);
	        	T = normalizeDataset(diag(s)*v');
    	    	[u,s,v] = svd(T,'econ');
        		s = diag(s);
        		displayLogItem(sprintf('\tCompleted re-normalizing dataset of subject %s!', subj), logFile);        	
	        end
        
    	    % Save out data to directory
        	saveOrigProjectedData(subj, v, s);
            
            if (useInitWarp == 1) && (continueFromWhereLeftOff==0 || checkWarpFileExists(subj, hems, 0)==0)
                % Then we also need to create the warped projected data for
                % each subject given the intial warps
                subWarpedTS = computeInterpOnSphere(diagMultPre(s, v'), coords.cartCoords, coords.neighbors, initWarp{subjNum});
                subWarpedTS = normalizeDataset(subWarpedTS);
                [u,s,v] = svd(subWarpedTS, 'econ');
                saveWarpedProjectedData(subj, v, diag(s));
                clear u s v subWarpedTS;
            end
	    end
	end
    displayLogItem('Completed pre-processing the datasets...', logFile);
    clear u s v p T;
%% END Data preprocessing
    

%% Determine the subject ordering
	if strcmp(firstPassAtlasMethod, 'build_up')
		% Then we need to determine the subject ordering
	    displayLogItem('Determining the alignment subject ordering...', logFile);
    	% Subjects will be ordered by decreasing average time-series correlation with the average of the remaining subjects
	    % Atlas will be leave-one-out atlas
	    atl = zeros(tLen, numNodes);
    	corrs = zeros(nsubj, 1);
    
	    % First, build up the atlas for the first subject
    	for subjNum = 2:nsubj
        	atl = atl + loadOrigTimeseries(subjects{subjNum}, startFromZScoredData);
	    end

    	for subjNum = 1:nsubj
        	subj = subjects{subjNum};

	        T = loadOrigTimeseries(subj, startFromZScoredData);
    	    if subjNum ~= 1
        	    % Then we need to subtract the subject from atl
            	atl = atl - T;
	        end

    	    for n = 1:numNodesPerHem
        	    magd = sqrt(sum(T(:, n).^2));
            	magd = magd*sqrt(sum(atl(:, n).^2));

	            if magd > 1e-8
    	            corrs(subjNum) = corrs(subjNum) + T(:, n)'*atl(:, n)/magd;
        	    end
	        end

    	    atl = atl + T;
	    end
    	% Now sort the correlations in descending order
	    [garbage, I] = sort(corrs, 'descend');
    	% Re-order the subjects based on correlation above
	    subjects = subjects(I);
        initWarp = initWarp(I);

        displayLogItem('Completed determining the alignment subject ordering', logFile);
    	clear atl garbage T corrs;
	elseif strcmp(firstPassAtlasMethod, 'average')
		% Then we don't need to determine the subject ordering
	else
		error(sprintf('Unknown value for firstPassAtlasMethod:  %s', firstPassAtlasMethod));
	end
    displayLogItem(['Subject ordering:  ', mat2str(cell2mat(subjects'))], logFile);
%% END Determine the subject ordering


%% Determine the starting point, if necessary
    startPass = 0;
    startSubjNum = 1;
	if continueFromWhereLeftOff == 1
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
            
    		zcDone = 1;
	    	if startSubjNum == 1 && startPass > 1
    			% Then we need to check if the zero correction has been done from the previous pass
    			for subjNum = 1:nsubj
    				if checkWarpFileExists(subjects{subjNum}, hems, startPass-1, 1) == 0
    					zcDone = 0;
	    				break;
    				end
    			end
				if zcDone == 1
					displayLogItem(['Zero correction was already completed for pass #', num2str(startPass-1)], logFile);
				else
					displayLogItem(['Zero correction did not complete from the previous pass.  Doing zero correction 	now...'], logFile);
					zeroCorrectWarps(startPass-1);
					displayLogItem(['Completed zero-correction from previous pass.'], logFile);
				end
	    	end
    
    		if zcDone == 1
	            % Check if warped projected data exists
    	        displayLogItem(['Checking if warped projected data exists...'], logFile);
				for subjNum = 1:nsubj
					subj = subjects{subjNum};
					dontWarp=0;
					if checkIfWarpedProjectedDataExists(subj) == 1
						displayLogItem(['Found warped projected data for ', subj, '!'], logFile);
						continue;
					else
						displayLogItem(['Need to recreate warped projected data for ', subj, '...'], logFile);
						if startSubjNum > subjNum
							displayLogItem(['Using warp from pass #', num2str(startPass), ' (not zero-corrected)'], logFile);
							warp = loadWarp(subj, startPass, 0);
						else
							if (startPass-1) > 0
								displayLogItem(['Using warp from pass #', num2str(startPass-1), ' (zero-corrected)'], logFile);
								warp = loadWarp(subj, startPass-1, 1);
							elseif startPass == 1
								displayLogItem(['Using warp from pass #0 (non zero-corrected)'], logFile);
								warp = loadWarp(subj, 0, 0);
							else
								% Then we are on pass #0 -- use original data -- don't do anything...
								dontWarp=1;
%								displayLogItem(['Using warp from pass #', num2str(startPass), ' (not zero-corrected)'], logFile);
%								warp = loadWarp(subj, startPass, 0);
							end
						end

						if dontWarp == 0
							% We need to regenerate the warped projected data for this subject...
							[subOrig.V, subOrig.sigma] = loadOrigProjectedData(subj);
							subWarpedTS = computeInterpOnSphere(diagMultPre(subOrig.sigma, transpose(subOrig.V)), coords.cartCoords, coords.neighbors, warp);
							subWarpedTS = normalizeDataset(subWarpedTS);
							[u,s,v] = svd(subWarpedTS, 'econ');
							subWarped.sigma = diag(s);
							subWarped.V = v;
							clear u s v subWarpedTS;
							saveWarpedProjectedData(subj, subWarped.V, subWarped.sigma);
							clear subWarped warp;
						end

						displayLogItem(['Completed recreating warped projected data for ', subj, '...'], logFile);
					end
				end
			end					
            displayLogItem(['Completed checking if warped projected data exists!'], logFile);
        end
    end
%% END Determine the starting point


%% 0th pass
    if startPass == 0
        displayLogItem('Running through pass #0...', logFile);

    	if strcmp(firstPassAtlasMethod, 'build_up')
	        if startSubjNum == 1
    	        % This is a special case -- we are just starting then
        	    % Save out the zero warp field for the first subject
            	saveWarp(zeros(3, numNodes), subjects{1}, 0, 0);
            
	            % Load in the data for the first subject   
                if useInitWarp == 1
                    subWarped = loadWarpedProjectedData(subjects{1});
                else
        	        subWarped = loadOrigProjectedData(subjects{1});
                end
        	    subWarped = normalizeDataset(transpose(subWarped));
            	[u,s,v] = svd(subWarped, 'econ');
	            subOrig.V = v;
    	        subOrig.sigma = diag(s);
        	    clear u s v subWarped;
            	% Save out "warped" data for subject #1
	            saveWarpedProjectedData(subjects{1}, subOrig.V, subOrig.sigma);
            
    	        % Now generate atlas
        	    atlas.V = subOrig.V;
            	atlas.sigma = subOrig.sigma;
	            clear subOrig;
    	    else
        	    % Then we need to create the atlas up to this point
            	% Load in all subjects that came before and average
	            % Also load in the warps to fill the subjWarps array
    	        displayLogItem('Re-creating atlas from previous run...', logFile);
        	    [atlas.V, atlas.sigma] = loadWarpedProjectedData(subjects{1});
            	for subjNum = 2:(startSubjNum-1)
                	[VA, sa] = loadWarpedProjectedData(subjects{subjNum});
	                [atlas.V, atlas.sigma, hasErrored] = updateGroupTemplate(atlas.V, atlas.sigma, (subjNum-1), VA, sa);
	                if hasErrored == 1
	                	whos
	                	error('Update group template didn''t have enough memory to complete');
	                end
    	        end
        	    clear VA sa;
            	displayLogItem('Completed re-creating atlas from previous run', logFile);
	        end

    	    if startSubjNum == 1
        	    startSubjNum = 2;
	        end
	    elseif strcmp(firstPassAtlasMethod, 'average')
	    	displayLogItem('Creating atlas...', logFile);
	    	% Then we need to build up the atlas to this point
		    atlasStarted = 0;
		    for subjNum = 1:(startSubjNum-1)
		    	if atlasStarted == 0
		    		[atlas.V, atlas.sigma] = loadWarpedProjectedData(subjects{subjNum});
		    		atlasStarted = 1;
		    	else
		    		% Then we need to add the subject into the current atlas
		    		[VA, sa] = loadWarpedProjectedData(subjects{subjNum});
		    		[atlas.V, atlas.sigma, hasErrored] = updateGroupTemplate(atlas.V, atlas.sigma, (subjNum-1), VA, sa);
	                if hasErrored == 1
	                	whos
	                	error('Update group template didn''t have enough memory to complete');
	                end
		    	end
            end
            % We want to skip startSubjNum  because we are building the
            % leave-one-out template for him!!!
			for subjNum = (startSubjNum+1):nsubj
                if useInitWarp == 1
                    [v,s] = loadWarpedProjectedData(subjects{subjNum});
                else
                	if startFromZScoredData == 1
                        [v,s] = loadOrigProjectedData(subjects{subjNum});
                    else
        				subOrig = normalizeDataset(transpose(loadOrigProjectedData(subjects{subjNum})));
            			[u,s,v] = svd(subOrig, 'econ');
                		s = diag(s);
                    end
                end
                
				if atlasStarted == 0
					atlas.V = v;
					atlas.sigma = s;
					atlasStarted = 1;
				else
					% Then we need to add the subject into the current atlas
					[atlas.V, atlas.sigma, hasErrored] = updateGroupTemplate(atlas.V, atlas.sigma, (subjNum-2), v, s);
	                if hasErrored == 1
	                	whos
	                	error('Update group template didn''t have enough memory to complete');
	                end
				end
			end
			clear u s v subOrig VA sa;
			displayLogItem('Completed creating atlas', logFile);
	    end

        % Now run through the 0th pass
        for subjNum = startSubjNum:nsubj
            subj = subjects{subjNum};

            % Load the data in
            [subOrig.V, subOrig.sigma] = loadOrigProjectedData(subj);   
            % Set the output file prefix
            alignLogFile = sprintf('%s/log_alignment_pass_0_subjNum_%d_subj_%s.txt', logDir, subjNum, subj);
            % Now run the algorithm
            displayLogItem(['Performing alignment of subject #', num2str(subjNum), ' (', subj, ') on pass #0'], logFile);
            displayLogItem(['A log of this alignment will be saved to ', alignLogFile], logFile);
            
            ncomponents = min(size(subOrig.V, 2), size(atlas.V, 2));
            atlasReg.V = atlas.V(:, 1:ncomponents);
            atlasReg.sigma = atlas.sigma(1:ncomponents);
            warp = connnorm_register(subOrig, atlasReg, coords, initWarp{subjNum}, regularization, maxResolution, alignLogFile);
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
            
          	if strcmp(firstPassAtlasMethod, 'build_up')
				if subjNum < nsubj
	                % Now update the atlas
    	            displayLogItem(['Updating group template by adding subject #', num2str(subjNum), ' into template'], logFile);
        	        [atlas.V, atlas.sigma, hasErrored] = updateGroupTemplate(atlas.V, atlas.sigma, (subjNum-1), subWarped.V, subWarped.sigma);
	                if hasErrored == 1
	                	whos
	                	error('Update group template didn''t have enough memory to complete');
	                end
        	        displayLogItem(['Completed updating group template!'], logFile);
	            else
	            	% Then we need to remove the first subject and add the last
    	        	displayLogItem(['Updating group template by adding subject #', num2str(subjNum-1), ' and removing subject #1'], logFile);
        	    	[subRemove.V, subRemove.sigma] = loadWarpedProjectedData(subjects{1});
            		[atlas.V, atlas.sigma, hasErrored] = updateGroupTemplate(atlas.V, atlas.sigma, (subjNum-1), subWarped.V, subWarped.sigma, subRemove.V, subRemove.sigma);
	                if hasErrored == 1
	                	whos
	                	error('Update group template didn''t have enough memory to complete');
	                end
            		clear subRemove;
            		displayLogItem(['Completed updating group template'], logFile);
            	end
			elseif strcmp(firstPassAtlasMethod, 'average')
				if subjNum < nsubj
        	    	% Now update the atlas
        	    	displayLogItem(['Updating group template by adding subject #', num2str(subjNum), ' and removing subject #', num2str(subjNum+1)], logFile);
        	    	if useInitWarp == 1
                        [subRemove.V, subRemove.sigma] = loadWarpedProjectedData(subjects{subjNum+1});
                    else
                        [subRemove.V, subRemove.sigma] = loadOrigProjectedData(subjects{subjNum+1});
                    end
        	    	[atlas.V, atlas.sigma, hasErrored] = updateGroupTemplate(atlas.V, atlas.sigma, (nsubj-1), subWarped.V, subWarped.sigma, subRemove.V, subRemove.sigma);
	                if hasErrored == 1
	                	whos
	                	error('Update group template didn''t have enough memory to complete');
	                end
        	    	clear subRemove;
        	    	displayLogItem(['Completed updating group template!'], logFile);        	    	
        	    else
        	    	displayLogItem(['Updating group template by adding subject #', num2str(subjNum), ' and removing subject #1'], logFile);
        	    	[subRemove.V, subRemove.sigma] = loadWarpedProjectedData(subjects{1});
        	    	[atlas.V, atlas.sigma, hasErrored] = updateGroupTemplate(atlas.V, atlas.sigma, (nsubj-1), subWarped.V, subWarped.sigma, subRemove.V, subRemove.sigma);
	                if hasErrored == 1
	                	whos
	                	error('Update group template didn''t have enough memory to complete');
	                end
        	    	clear subRemove;
        	    	displayLogItem(['Completed updating group template!'], logFile);
        	    end
        	else
        	    error('Unknown firstPassAtlasMethod: %s', firstPassAtlasMethod);
        	end
            clear subWarped;
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
        atlasStarted = 0;
        subjCount = 0;
        for subjNum = 1:nsubj
        	if subjNum == startSubjNum; continue; end;
        	
        	if atlasStarted == 0
        		[atlas.V, atlas.sigma] = loadWarpedProjectedData(subjects{subjNum});
        		atlasStarted = 1;
        	else
        		[VA, sa] = loadWarpedProjectedData(subjects{subjNum});
        		[atlas.V, atlas.sigma, hasErrored] = updateGroupTemplate(atlas.V, atlas.sigma, subjCount, VA, sa);
	                if hasErrored == 1
	                	whos
	                	error('Update group template didn''t have enough memory to complete');
	                end
        	end
        	
        	subjCount = subjCount + 1;
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
            alignLogFile = sprintf('%s/log_alignment_pass_%d_subjNum_%d_subj_%s.txt', logDir, passNum, subjNum, subj);
            % Load the current warp file 
            if passNum > 1
                warp = loadWarp(subj, passNum-1, 1);
            else
                warp = loadWarp(subj, passNum-1, 0);
            end
            % Now run the algorithm            
            displayLogItem(['Performing alignment of subject #', num2str(subjNum), ' (', subj, ') on pass #', num2str(passNum)], logFile);
            displayLogItem(['A log of this alignment will be saved to ', alignLogFile], logFile);
            ncomponents = size(subOrig.V, 2);
            atlasReg.V = atlas.V(:, 1:ncomponents);
            atlasReg.sigma = atlas.sigma(1:ncomponents);
            warp = connnorm_register(subOrig, atlasReg, coords, warp, regularization, maxResolution, alignLogFile);
            displayLogItem(['Completed alignment of subject #', num2str(subjNum), ' (', subj, ') on pass #', num2str(passNum)], logFile);
            
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
            	[atlas.V, atlas.sigma, hasErrored] = updateGroupTemplate(atlas.V, atlas.sigma, (nsubj-1), subWarped.V, subWarped.sigma, subNext.V, subNext.sigma);
	                if hasErrored == 1
	                	whos
	                	error('Update group template didn''t have enough memory to complete');
	                end
            	clear subNext;
            end
        end

		displayLogItem('Zero correcting the warps, and re-adjusting the atlas...', logFile);
		zeroCorrectWarps(passNum);
		displayLogItem('Completed zero correcting the warps, and re-adjusting the atlas', logFile);
        
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
    





%% FUNCTION zeroCorrectWarps
%	This function zero corrects warps and resets atlas for 1st subject
	function zeroCorrectWarps(passNum)
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
					[atlas.V, atlas.sigma, hasErrored] = updateGroupTemplate(atlas.V, atlas.sigma, (subjNum - 2), subWarped.V, subWarped.sigma);
	                if hasErrored == 1
	                	whos
	                	error('Update group template didn''t have enough memory to complete');
	                end
				end
			end                
        end        
	end
%% END FUNCTION zeroCorrectWarps


%% FUNCTION getOrigTimeSeriesFilename
    function fName = getOrigTimeSeriesFilename(subj, hem)
        fName = sprintf('%s/%s_%s_%dmm_fwhm%d_%s_on_%s.niml.dset', outDir, subj, hem, mmRes, fwhm, experiment, fs_surf);
    end
%% END FUNCTION getOrigTimeSeriesFilename


%% FUNCTION getOrigProjectedDataFilenames
    function [V_fName, s_fName] = getOrigProjectedDataFilenames(subj)
        V_fName = sprintf('%s/V_%s_orig.bin', tmpDir2, subj);
        s_fName = sprintf('%s/s_%s_orig.bin', tmpDir2, subj);
    end
%% END FUNCTION getOrigProjectedDataFilenames


%% FUNCTION checkIfOrigProjectedDataExists
	function e = checkIfOrigProjectedDataExists(subj)
		[V_fName, s_fName] = getOrigProjectedDataFilenames(subj);
		if exist(V_fName) == 2 && exist(s_fName) == 2
			e = 1;
			return;
		end
		e = 0;
	end
%% END FUNCTION checkIfOrigProjectedDataExists


%% FUNCTION getWarpedProjectedDataFilenames
    function [V_fName, s_fName] = getWarpedProjectedDataFilenames(subj)
        V_fName = sprintf('%s/V_%s_warped.bin', tmpDir2, subj);
        s_fName = sprintf('%s/s_%s_warped.bin', tmpDir2, subj);
    end
%% END FUNCTION getTmpTimeseriesFilenames

%% FUNCTION checkIfWarpedProjectedDataExists
	function e = checkIfWarpedProjectedDataExists(subj)
		[V_fName, s_fName] = getWarpedProjectedDataFilenames(subj);
		if exist(V_fName) == 2 && exist(s_fName) == 2
			e = 1;
			return;
		end
		e = 0;
	end
%% END FUNCTION checkIfWarpedProjectedDataExists



%% FUNCTION loadOrigTimeseries
    function T = loadOrigTimeseries(subj, zscore)
    	if ~exist('zscore') zscore = 0; end;
    	
        T = zeros(tLenFull, nhem*numNodesPerHem, opsDataType);

        n = 1;
        for hemNum = 1:nhem
            hem = hems{hemNum};
            T(:, n:(n+numNodesPerHem-1)) = read_niml_binary(getOrigTimeSeriesFilename(subj, hem), opsDataType, saveDataType);

            n = n + numNodesPerHem;
        end
        
		if useTWin == 1
			T = T(tWin, :);
		end
		
		T = centerDataset(T);
		if zscore == 1
			T = normalizeDataset(T);
		end
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
        try
	        numEigs = (d.bytes)/saveDataBytes;
    	catch
    		error(sprintf('s_fName directory:  %s', s_fName));
    	end
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
    function e = checkWarpFileExists(subj, hems, passNum, zeroCorrected)
    	if nargin == 3; zeroCorrected = 0; end;
    	
        e = 1;
        for hemNum = 1:nhem
            hem = hems{hemNum};

            ex = exist(getWarpFilename(subj, hem, passNum, zeroCorrected));
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

        warpFilename = sprintf('%s/warp_%s_%s_%d%s.bin', tmpDir, subj, hem, passNum, corrStr);
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
            
            fp = fopen(sprintf('%s/warp_%s_%s_final.bin', hemDir, subj, hem), 'wb');
            wh = warp(:, n:(n+numNodesPerHem-1));
            fwrite(fp, wh, saveCoordsType);
            fclose(fp);
            
            fp = fopen(sprintf('%s/standard%dmm_%s_%s.connnorm.asc', hemDir, mmRes, subj, hem), 'w');
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
        	subj = subjects{subjNum};
            for hemNum = 1:nhem
                hem = hems{hemNum};
                hemDir = hemOutDirs{hemNum};
                
                ex = exist(sprintf('%s/warp_%s_%s_final.bin', hemDir, subj, hem));
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
