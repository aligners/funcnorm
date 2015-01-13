function warps = script_align_subject_to_averageTimeSlice(tLen, lambda, pctEnergy)

subjects = {'cb', 'dm', 'hj', 'kd', 'kl', 'mh', 'ph', 'rb', 'se', 'sm'};
hems = {'lh', 'rh'};
experiment = 'raidersP1';
mmRes = 2;
fwhm = 0;
fs_surf = 'sphere.reg';
outDir = '../surf_data/connnorm_out_data/';
lambda_areal = lambda;
lambda_metric = lambda;
maxResolution = 3;
numPasses = 3;
logFile = 'connnorm_log.log';

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
    
  
	alignDir = [outDir, 'alignmentSingleSubjToAvgTimeSlice_', experiment, '_', num2str(lambda_areal), '_', num2str(pctEnergy), '_tLen', num2str(tLen), '/'];
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
    	displayLogItem(sprintf('\t%s', ['Pct Energy:  ', num2str(pctEnergy)]), logFile);
    	displayLogItem(sprintf('\t%s', ['Max resolution:  ', num2str(maxResolution)]), logFile);
    	displayLogItem(sprintf('\t%s', ['Num passes:  ', num2str(numPasses)]), logFile);
	
    

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
    TRs = 1:tLen;
%    tLen = getNumTimepoints(getOrigTimeSeriesFilename(subjects{1}, hems{1}));
%	if partNum == 1
%		TRs = 1:floor(tLen/2);
%		tLen = length(TRs);
%	else
%		TRs = (floor(tLen/2)+1):tLen;
%		tLen = length(TRs);
%	end
	displayLogItem(['Number of nodes per hemisphere:  ', num2str(numNodesPerHem)], logFile);
	displayLogItem(['Number of triangles per hemisphere:  ', num2str(numTrisPerHem)], logFile);
	displayLogItem(['Number of timepoints:  ', num2str(tLen)], logFile);
	displayLogItem(['Using TR #', num2str(TRs(1)), ' to #', num2str(TRs(end))], logFile);
	displayLogItem(['Number of hemispheres:  ', num2str(nhem)], logFile);
	
%% Setup regularization
    regularization.mode = 'metric_and_areal_distortion';
    displayLogItem('Using metric and areal distortion as warp regularization', logFile);
    displayLogItem(sprintf('\tLambda areal:  %f', lambda_areal), logFile);
    displayLogItem(sprintf('\tLambda metric:  %f', lambda_metric), logFile);
    regularization.lambda_areal = numNodes*lambda_areal;
    regularization.lambda_metric = numNodes*lambda_metric;
%% END Setup regularization


%% Data preprocessing
    displayLogItem('Pre-processing the datasets...', logFile);      
%    displayLogItem('WARNING:  TEMPORARILY DISABLED!!!', logFile);  
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

		subOrigTS = normalizeDataset(diag(s)*v');
		[u,s,v] = svd(subOrigTS, 'econ');
		s = diag(s);
		clear u;
		
        % Save out data to directory
        saveOrigProjectedData(subj, v, s);
    end
    displayLogItem('Completed pre-processing the datasets...', logFile);
    clear s v p T;
%% END Data preprocessing
    


for subjNum = 1:nsubj
	%% Generate the atlas at this point
	displayLogItem('Generating the atlas...', logFile);
	atlCount = 0;
	for atlSubjNum = 1:nsubj
		if atlSubjNum == subjNum
			continue;
		end

		atlSubj = subjects{atlSubjNum};	
		[atlSubOrig.V, atlSubOrig.sigma] = loadOrigProjectedData(atlSubj);

		if atlCount == 0
			atlas.V = atlSubOrig.V;
			atlas.sigma = atlSubOrig.sigma;
		else
			[atlas.V, atlas.sigma] = updateGroupTemplate(atlas.V, atlas.sigma, atlCount, atlSubOrig.V, atlSubOrig.sigma);
		end
		atlCount = atlCount + 1;
	end
	displayLogItem('Completed generating the atlas!', logFile);

	subj = subjects{subjNum};
	[subOrig.V, subOrig.sigma] = loadOrigProjectedData(subj);

	% Perform the alignment
	alignLogFile = [logDir, 'log_alignment_subjNum_', num2str(subjNum), '_subj_', subj, '.txt'];
	% Now run the algorithm
	displayLogItem(['Performing alignment of subject #', num2str(subjNum), ' (', subj, ')'], logFile);
	displayLogItem(['A log of this alignment will be saved to ', alignLogFile], logFile);
	ncomponents = min(size(subOrig.V, 2), size(atlas.V, 2));
	atlasReg.V = atlas.V(:, 1:ncomponents);
	atlasReg.sigma = atlas.sigma(1:ncomponents);
	warp = connnorm_register(subOrig, atlasReg, coords, zeros(3, numNodes, opsCoordsType), regularization, maxResolution, alignLogFile);
	displayLogItem(['Completed alignment of subject #', num2str(subjNum), ' (', subj, ')'], logFile);

	% Now we need to output the warp file
	saveWarp(warp, subj);
	
%	displayLogItem(['Updating the template...'], logFile);
%	% Update the template
%	if subjNum < nsubj
%		[subNext.V, subNext.sigma] = loadOrigProjectedData(subjects{subjNum+1});
%		[atlas.V, atlas.sigma] = updateGroupTemplate(atlas.V, atlas.sigma, (nsubj-1), subOrig.V, subOrig.sigma, subNext.V, subNext.sigma);
%		clear subNext;
%	end
%	displayLogItem(['Completed updating the template!'], logFile);
end
    








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
        T = cast(T(TRs, :), opsDataType);
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