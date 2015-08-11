function funcnorm_prepare_wrapper(subjects, hems, experiment, volDataDir, suffix, ...
    extensions, inpSubjDir, saveSubjDir, preprocess_movie, mmRes, fwhm, ...
    fs_surf, outDir, logFile, warp_experiment, alignDir, inpUnixPath, functype)

% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

	if ~exist('alignDir') || isempty(alignDir)
		specifyAlignDir = 0;
	else
		specifyAlignDir = 1;
	end
	
    if ~exist('logFile') || isempty(logFile)
        % Don't append to a logfile
        logFile = '';
    end
    
	if ~exist('warp_experiment') || isempty(warp_experiment)
		warpExperimentStr = '';
	else
 		warpExperimentStr = [' -warp_experiment ',warp_experiment];
    end
   	
    if ~exist('inpUnixPath') || isempty(inpUnixPath)
        inpUnixPath=getenv('PATH');
    end
    
    if ~exist('functype') || isempty(functype)
        functype='connnorm'
    end
    
    nsubj = length(subjects);
    nhem = length(hems);
            
	displayLogItem('Preparing the data with funcnorm_prepare.sh script...', logFile);

    % Setup environment funcnorm_prepare.sh expects
    % FREESURFER_HOME
    [retval fs_bin]=unix('which mris_convert');
    if retval
        error('mris_convert not found');
    end
    fs_bin=fs_bin(1:length(fs_bin)-1); % chop off end of line
    fs_dir=[fs_bin,'/../..'];
    setenv('FREESURFER_HOME', fs_dir);
    % unset LD_LIBRARY_PATH, which messes up FreeSurfer on Linux
    c=computer;
    islinux=strcmp(c(1:4),'GLNX');
    if islinux
        orig_ldpath=getenv('LD_LIBRARY_PATH');
        setenv('LD_LIBRARY_PATH');
    end
    % SUBJECTS_DIR
    setenv('SUBJECTS_DIR', inpSubjDir);
    % PATH
    setenv('PATH',inpUnixPath);
    
    if preprocess_movie
        preprocessFlagStr = ' -preprocess_movie';
    else
        preprocessFlagStr = '';
    end
    
    if specifyAlignDir
    	alignDirStr = [' -alignDir ', alignDir];
    else
    	alignDirStr = '';
    end

    for subjNum = 1:nsubj
        subj = subjects{subjNum};

        for hemNum = 1:nhem
            hem = hems{hemNum};

            displayLogItem(['Preparing ', subj, ' (', hem, ') with funcnorm_prepare.sh'], logFile);
            unixStr = ['sh ./funcnorm_prepare.sh ',...
                            '-sub ', subj,...
                            ' -hem ', hem,...
                            ' -inputdir ', volDataDir, ...
                            ' -experiment ', experiment, ...
                            warpExperimentStr, ...
                            ' -suffix ', suffix, ...
                            ' -subjects_dir ', saveSubjDir(1:(end-1)),...
                            ' -outputdir ', outDir(1:(end-1)),...
                            ' -fs_surf ', fs_surf,...
                            ' -mm ', num2str(mmRes),...
                            ' -fwhm ', num2str(fwhm),...
                            alignDirStr,...
                            preprocessFlagStr,...
                            ' -verbosity 2',...
                            ' -functype ' functype
                      ];                  
                  
            disp(unixStr);
            retval = unix(unixStr);
            if retval
                error('problem running funcnorm_prepare');
            end
        end
    end
    % put back LD_LIBRARY_PATH
    if islinux
        setenv('LD_LIBRARY_PATH',orig_ldpath)
    end
    
    displayLogItem('Completed preparing the data with funcnorm_prepare.sh script', logFile);
