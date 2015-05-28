function warps = funcnorm(subjects, hems, varargin)
% FUNCTION warps = funcnorm(subjects, hems, varargin)
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
% FUNCTION warps = funcnorm(subjects, hems, experiment, mmRes, fwhm, fs_surf, outDir, ...
%   lambda_areal, lambda_metric, maxResolution, numPasses [, logFile])
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
%   maxResolution:  multiplier specifying the maximum resolution for multi-resolution alignment
%       The alignment begins at a resolution of maxResolution*mmRes, where
%       mmRes is specified above as the initial cortical mesh resolution
%   numPasses:  the number of alignment passes through the subjects
%   logFile (optional):  filename to store log items (relative to outDir)
%
%
% FUNCTION warps = funcnorm(subjects, hems, experiment, volDataDir, ...
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
%   maxResolution:  multiplier specifying the maximum resolution for multi-resolution alignment
%       The alignment begins at a resolution of maxResolution*mmRes, where
%       mmRes is specified above as the initial cortical mesh resolution
%   numPasses:  the number of alignment passes through the subjects
%   logFile (optional):  filename to store log items (relative to outDir)

% Copyright (c) 2008 Bryan Conroy, Jim Haxby, Peter Ramadge, Mert Rory Sabuncu, and Ben Singer


    if nargin == 11 || nargin == 12
        % FUNCTION warps = funcnorm(subjects, hems, experiment, mmRes, fwhm, fs_surf, outDir, ...
        %   lambda_areal, lambda_metric, maxResolution, numPasses [, logFile])
        
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
        numPasses = varargin{9};
        
        if nargin == 11
            useLog = 0;
        else
            useLog = 1;
            logFile = varargin{10};
        end
    elseif nargin == 17 || nargin == 18
        % FUNCTION warps = funcnorm(subjects, hems, experiment, volDataDir, ...
        %   suffix, extensions, inpSubjDir, saveSubjDir, preprocess_movie, ...
        %   mmRes, fwhm, fs_surf, outDir, lambda_areal, lambda_metric, maxResolution, ... 
        %   numPasses [, logFile])
        
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
        maxResolution = varargin{14};
        numPasses = varargin{15};
        
        if nargin == 17
            useLog = 0;
        else
            useLog = 1;
            logFile = varargin{16};
        end
    else
        error('Improper call to funcnorm.m  Type ''help funcnorm'' for details.');
    end

    saveDataBytes = 4; % Number of bytes of datatype to store data (e.g. float32 = 4 bytes)
    saveDataType = 'float32';
    opsDataBytes = 4; % Number of bytes to perform all operations in (e.g. single precision = 4 bytes)
    opsDataType = 'single';

    if outDir(end) ~= '/'
        outDir = [outDir, '/'];
    end
    
    alignDir = [outDir, 'alignment_', experiment, '/'];
    if isdir(alignDir) == 0
        mkdir(alignDir);
    end
    
    logDir = [alignDir, 'log/'];
        if isdir(logDir) == 0
            mkdir(logDir);
        end
    warpsDir = [alignDir, 'warps/'];
        if isdir(warpsDir) == 0
            mkdir(warpsDir);
        end
    if useLog == 1
        logFile = [logDir, logFile];
    else
        logFile = '';
    end
    
    nsubj = length(subjects); % Number of subjects
    nhem = length(hems); % Number of hemispheres

    displayLogItem('Beginning funcnorm algorithm...');
    
    
%% Prepare the data if necessary
    if runFuncnormPrepare == 1
        funcnorm_prepare_wrapper(subjects, hems, experiment, volDataDir, suffix, ...
            extensions, inpSubjDir, saveSubjDir, preprocess_movie, mmRes, fwhm, ...
            fs_surf, outDir, logFile)
    end
%% END Prepare the data if necessary
    

%% Parse the surface file and setup coordinates structure
    displayLogItem('Parsing the surface file...', logFile);
    surfFile = [outDir, 'standard', num2str(mmRes), 'mm_', fs_surf, '.asc'];
    [numNodesPerHem, coords] = parseSurfaceFile(surfFile, opsDataType);
    numNodes = numNodesPerHem*nhem;
    numTrisPerHem = size(coords.triangles, 2);
    numTris = nhem*numTrisPerHem;
    for j = 1:numNodesPerHem
        coords.cartCoords(:, j) = coords.cartCoords(:, j)/sqrt(sum(coords.cartCoords(:, j).^2));
    end
    % Replicate everything for all hemispheres
    if nhem > 1
        coords = generateMultiHemisphereCoordsDataset(coords, nhem);
    end
    displayLogItem('Completed parsing the surface file', logFile);
%% END Parse the surface file and setup coordinates structure

    % Determine the number of time-points
    tLen = getNumTimepoints(getOrigTimeSeriesFilename(subjects{1}, hems{1}));

%% Determine the subject ordering
    displayLogItem('Determining the alignment subject ordering...', logFile);
    % First subject is the one with highest average inter-subject correlation
    atl = zeros(tLen, numNodes);
    corrs = zeros(nsubj, 1);
    
    % First, build up the atlas
    for subjNum = 2:nsubj
        atl = atl + loadOrigTimeseries(subjects{subjNum});
    end

    for subjNum = 1:nsubj
        subj = subjects{subjNum};

        T = loadOrigTimeseries(subj);
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
    [garbage, I] = sort(corrs, 'descend');
    % Re-order the subjects
    subjects = subjects(I);
    displayLogItem('Completed determining the alignment subject ordering', logFile);
    displayLogItem(['Subject ordering:  ', mat2str(cell2mat(subjects'))], logFile);
    clear atl;
%% END Determine the subject ordering

%% Setup regularization
    regularization.mode = 'metric_and_areal_distortion';
    regularization.lambda_areal = numNodes*lambda_areal;
    regularization.lambda_metric = numNodes*lambda_metric;
%% END Setup regularization

%% Setup output directory structure
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

%% Determine the starting point
    startPass = 0;
    startSubjNum = 1;
    if determineStartPoint == 1
        % If the directory already exists, then try to determine the starting
        % point
        if checkWarpCompleted == 1
            displayLogItem('This atlas construction has completed already.');
            return
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
%% END Determine the starting point

%% 0th pass
    if startPass == 0
        displayLogItem('Running through pass #0...', logFile);

        if startSubjNum == 1
            % This is a special case -- we are just starting then
            % Load in the data for the first subject and label it as atlas
            % Save out the zero warp field for the first subject
            saveWarp(zeros(3, numNodes), subjects{1}, 0, 0);
            atlasTS = loadOrigTimeseries(subj);
            % Save out "warped" timeseries for subject #1
            saveTmpTimeseries(atlasTS, subjects{1});
        else
            % Then we need to create the atlas up to this point
            % Load in all subjects that came before and average
            % Also load in the warps to fill the subjWarps array
            displayLogItem('Re-creating atlas from previous run...', logFile);
            atlasTS = zeros(tLen, numNodes);
            for j = 1:(startSubjNum-1)
                atlasTS = atlasTS + loadTmpTimeseries(subjects{j});
            end
            displayLogItem('Completed re-creating atlas from previous run', logFile);
        end

        if startSubjNum == 1
            startSubjNum = 2;
        end

        % Now run through the 0th pass
        for subjNum = startSubjNum:nsubj
            subj = subjects{subjNum};

            % Load the data in
            TS = loadOrigTimeseries(subj);   
            % Set the output file prefix
            alignLogFile = [logDir, 'log_alignment_pass_0_subjNum_', num2str(subjNum), '_subj_', subj, '.txt'];
            % Now run the algorithm
            displayLogItem(['Performing alignment of subject #', num2str(subjNum), ' (', subj, ') on pass #0'], logFile);
            displayLogItem(['A log of this alignment will be saved to ', alignLogFile], logFile);
            
            warp = funcnorm_register(atlasTS, TS, coords, zeros(3, numNodes, opsDataType), regularization, maxResolution, alignLogFile);
            displayLogItem(['Completed alignment of subject #', num2str(subjNum), ' (', subj, ') on pass #0'], logFile);

            % Now we need to output the warp file
            saveWarp(warp, subj, 0, 0);
            % Now interpolate the data using the warp
            TS = computeInterpOnSphere(TS, coords.cartCoords, coords.neighbors, warp);
            % Write out the temp dataset to file
            saveTmpTimeseries(TS, subj);
            % Now average with the atlas
            atlasTS = TS + atlasTS;
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
        atlasTS = zeros(tLen, numNodes);
        for subjNum = 1:nsubj
            subj = subjects{subjNum};

            atlasTS = atlasTS + loadTmpTimeseries(subj);
        end
        displayLogItem('Completed re-creating atlas from previous run', logFile);
    end
%% END 0th pass

%% Remaining passes
    % Now run through the remaining passes
    for passNum = startPass:numPasses
        displayLogItem(['Running through pass #', num2str(passNum), '...'], logFile);
        for subjNum = startSubjNum:nsubj
            subj = subjects{subjNum};

            % Load in the (warped) time series
            TS = loadTmpTimeseries(subj);
            % Subtract it from the atlas
            atlasTS = atlasTS - TS;
            clear TS;

            % Load in the original timeseries
            TS = loadOrigTimeseries(subj);
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
            warp = funcnorm_register(TS, atlasTS, coords, warp, regularization, maxResolution, alignLogFile);
            displayLogItem(['Completed alignment of subject #', num2str(subjNum), ' (', subj, ') on pass #', num2str(passNum)], logFile);
            
            % Now we need to output the warp file
            saveWarp(warp, subj, passNum, 0);
            % Now interpolate the data using the warp
            TS = computeInterpOnSphere(TS, coords.cartCoords, coords.neighbors, warp);
            % Write out the temp dataset to file
            saveTmpTimeseries(TS, subj);
            % Now average with the atlas
            atlasTS = TS + atlasTS;    
        end

        if passNum > 0
            displayLogItem('Zero correcting the warps...', logFile);
            % Then zero correct the warps
            for subjNum = 1:nsubj
                subj = subjects{subjNum};
                warps{subjNum} = loadWarp(subj, passNum, 0);
            end

            warps_zero = computeZeroCorrection(coords.cartCoords, warps);
            for subjNum = 1:nsubj
                subj = subjects{subjNum};

                saveWarp(warps_zero{subjNum}, subj, passNum, 1);
            end        
            displayLogItem('Completed zero correcting the warps', logFile);
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


%% FUNCTION getTmpTimeseriesFilename
    function fName = getTmpTimeseriesFilename(subj, hem)
        fName = [tmpDir, subj, '_', hem, '.bin'];
    end
%% END FUNCTION getTmpTimeseriesFilename


%% FUNCTION loadOrigTimeseries
    function T = loadOrigTimeseries(subj)
        T = zeros(tLen, nhem*numNodesPerHem, opsDataType);

        n = 1;
        for hemNum = 1:nhem
            hem = hems{hemNum};
            T(:, n:(n+numNodesPerHem-1)) = read_niml_binary(getOrigTimeSeriesFilename(subj, hem));

            n = n + numNodesPerHem;
        end
    end
%% END FUNCTION loadOrigTimeseries


%% FUNCTION loadTmpTimeseries
    function T = loadTmpTimeseries(subj)
        T = zeros(tLen, nhem*numNodesPerHem, opsDataType);

        n = 1;
        for hemNum = 1:nhem
            hem = hems{hemNum};

            fp = fopen(getTmpTimeseriesFilename(subj, hem), 'rb');
            T(:, n:(n+numNodesPerHem-1)) = fread(fp, [tLen, numNodesPerHem], saveDataType);
            fclose(fp);

            n = n + numNodesPerHem;
        end
    end
%% END FUNCTION loadTmpTimeseries

    
%% FUNCTION saveTmpTimeseries
    function saveTmpTimeseries(T, subj)
        n = 1;
        for hemNum = 1:nhem
            hem = hems{hemNum};

            fp = fopen(getTmpTimeseriesFilename(subj, hem), 'wb');
            fwrite(fp, T(:, n:(n+numNodesPerHem-1)), saveDataType);
            fclose(fp);

            n = n + numNodesPerHem;
        end
    end
%% END FUNCTIOn saveTmpTimeseries


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
        warp = zeros(3, nhem*numNodesPerHem, opsDataType);

        n = 1;
        for hemNum = 1:nhem
            hem = hems{hemNum};

            fp = fopen(getWarpFilename(subj, hem, passNum, zeroCorrected), 'rb');
            warp(:, n:(n+numNodesPerHem-1)) = fread(fp, [3, numNodesPerHem], saveDataType);
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
            fwrite(fp, warp(:, n:(n+numNodesPerHem-1)), saveDataType);
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
            fwrite(fp, wh, saveDataType);
            fclose(fp);
            
            fp = fopen([hemDir, 'standard', num2str(mmRes), 'mm_', subj, '_', hem, '.funcnorm.asc'], 'w');
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
%% END funcnorm