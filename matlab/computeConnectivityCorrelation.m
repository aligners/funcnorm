function computeInterSubjectCorrelationResults(subjects, hems, outDir, train_experiment, test_experiment, mmRes, fwhm, fs_surf)
% FUNCTION computeInterSubjectCorrelationResults(subjects, hems, outDir, ...
%   train_experiment, test_experiment, mmRes, fwhm, fs_surf)
%
% For anatomical results, set train_experiment = 'anat'

sep = sprintf('\t');
logFile = '';
incTime = 0;


if outDir(end) ~= '/'
    outDir = [outDir, '/'];
end

resultsDir = [outDir, 'alignment_', train_experiment, '/results/'];
if isdir(resultsDir) == 0
    mkdir(resultsDir);
end

nsubj = length(subjects);
nhem = length(hems);

tLen = getNumTimepoints(getOrigTimeSeriesFilename(subjects{1}, hems{1}, test_experiment));
numNodesPerHem = getNumNodes(getOrigTimeSeriesFilename(subjects{1}, hems{1}, test_experiment));
numNodes = nhem*numNodesPerHem;


corrs = zeros(nsubj, numNodes);
loopStep = 2000;
numLoops = ceil(numNodes/loopStep);


for k = 1:loopStep:numNodes
	displayLogItem(['Working on loop #', num2str((k-1)/loopStep + 1), ' out of ', num2str(numLoops), '...'], logFile, incTime);
	
	kMax = min(k+loopStep-1, numNodes);
	kInd = k:kMax;
	currNumNodes = length(kInd);
	
	atl = zeros(numNodes, currNumNodes);

    for subjNum = 1:nsubj
        c{subjNum} = zeros(numNodes, currNumNodes);
    end

	displayLogItem([sep, 'Computing leave-one-out-atlas...'], logFile, incTime);
    for subjNum = 1:nsubj
        subj = subjects{subjNum};
        
        if strcmp(train_experiment, 'anat')
        	T = loadOrigTimeseries(subj, test_experiment);
        else
        	T = loadWarpedTimeseries(subj, test_experiment, train_experiment);
        end

		% Now center the datasets
		T = normalizeDataset(T);
		
		c{subjNum} = T'*T(:, kInd);

		if subjNum > 1
			atl = atl + c{subjNum}/(nsubj-1);
		end
    end
	displayLogItem([sep, 'Completed computing leave-one-out-atlas!'], logFile, incTime);

	displayLogItem([sep, 'Computing correlations...'], logFile, incTime);
	for subjNum = 1:nsubj
		subj = subjects{subjNum};

		if subjNum > 1
			% Remove the subject from the atlas, and add in previous subject
			atl = atl + (c{subjNum-1} - c{subjNum})/(nsubj - 1);
		end
		
		sC = c{subjNum};
		for n = 1:currNumNodes
			R = corrcoef(sC(:, n), atl(:, n));
			corrs(subjNum, kInd(n)) = R(2);
		end
	end
	displayLogItem([sep, 'Completed computing correlations!'], logFile, incTime);

    displayLogItem(['Completed loop #', num2str((k-1)/loopStep + 1), ' out of ', num2str(numLoops), '!'], logFile, incTime);
end
% corrcoef stupidly sets to NaN if vector has zero norm
corrs(find(isnan(corrs))) = 0;

% Write out the connectivity correlation results for each subject
displayLogItem('Writing out connectivity correlation results for each subject...', logFile, incTime);
for subjNum = 1:nsubj
	subj = subjects{subjNum};
	n = 1;
	for hemNum = 1:nhem
		hem = hems{hemNum};
		
		outStr = [resultsDir, train_experiment, 'Align_connectivitycorrelation_', test_experiment, '_', subj, '_', hem, '.dset'];
		displayLogItem([sep, 'Connectivity correlation results for subject (', subj, '), ', hem, ' hemisphere saved to:  ', outStr], logFile, incTime);
		fp = fopen(outStr, 'w');
		fprintf(fp, '%f\n', corrs(subjNum, n:(n+numNodesPerHem-1)));
		fclose(fp);
		
		n = n + numNodesPerHem;
	end
end
displayLogItem('Completed writing out connectivity correlation results for each subject!', logFile, incTime);
	
% Get the average connectivity correlation at each node
displayLogItem('Writing out average connectivity correlation results for each hemisphere...', logFile, incTime);
corrs_avg = sum(corrs)/nsubj;
% Get the total average connectivity correlation of each hemisphere
totalAvgCorrPerHem = zeros(1, nhem);
n = 1;
for hemNum = 1:nhem
    hem = hems{hemNum};
    totalAvgCorrPerHem(hemNum) = sum(corrs_avg(n:(n+numNodesPerHem-1)))/numNodesPerHem;
    
    outStr = [resultsDir, train_experiment, 'Align_avg_connectivitycorrelation_', test_experiment, '_', hem, '.dset'];
    displayLogItem([sep, 'Average connectivity correlation results for ', hem, ' hemisphere saved to:  ', outStr], logFile, incTime);
    fp = fopen(outStr, 'w');
    fprintf(fp, '%f\n', corrs_avg(n:(n+numNodesPerHem-1)));
    fclose(fp);
    
    n = n + numNodesPerHem;
end
displayLogItem('Completed writing out average connectivity correlation results for each hemisphere!', logFile, incTime);



% Get the total average inter-subject correlation across all hemispheres
%anat_totalAvgCorr = sum(anat_totalAvgCorrPerHem)*numNodesPerHem/numNodes;
%displayLogItem(sprintf('\n'), logFile, incTime);
%displayLogItem(['Anatomical Alignment Results (on ', test_experiment, '):'], logFile, incTime);
%    displayLogItem([sep, 'Average inter-subject correlation (per hemisphere):  '], logFile, incTime);
%        for hemNum = 1:nhem
%            hem = hems{hemNum};
%            displayLogItem([sep, sep, hem, ':  ', num2str(anat_totalAvgCorrPerHem(hemNum))], logFile, incTime);
%        end
%    displayLogItem([sep, 'Total average inter-subject correlation:  ', num2str(anat_totalAvgCorr)], logFile, incTime);

%displayLogItem(sprintf('\n'), logFile, incTime);
%displayLogItem(['Functional Alignment Results (on ', test_experiment, '):'], logFile, incTime);
%    displayLogItem([sep, 'Average inter-subject correlation (per hemisphere):  '], logFile, incTime);
%        for hemNum = 1:nhem
%            hem = hems{hemNum};
%            displayLogItem([sep, sep, hem, ':  ', num2str(func_totalAvgCorrPerHem(hemNum))], logFile, incTime);
%        end
%    displayLogItem([sep, 'Total average inter-subject correlation:  ', num2str(func_totalAvgCorr)], logFile, incTime);





    
%% FUNCTION getOrigTimeSeriesFilename
function fName = getOrigTimeSeriesFilename(subj, hem, experiment)
    fName = [outDir, subj, '_', hem, '_', num2str(mmRes), 'mm_fwhm', num2str(fwhm), '_', experiment, '_on_', fs_surf, '.niml.dset'];
end
%% END FUNCTION getOrigTimeSeriesFilename


%% FUNCTION getWarpedTimeseriesFilename
function fName = getWarpedTimeSeriesFilename(subj, hem, test_experiment, train_experiment)
    fName = [outDir, subj, '_', hem, '_', num2str(mmRes), 'mm_fwhm', num2str(fwhm), '_', test_experiment, '_warpedby_', train_experiment, '_on_', fs_surf, '.niml.dset'];
end
%% END FUNCTION getWarpedTimeseriesFilename


%% FUNCTION loadOrigTimeseries
function T = loadOrigTimeseries(subj, experiment)
    T = zeros(tLen, numNodes);

    n = 1;
    for hemNum = 1:nhem
        hem = hems{hemNum};
        T(:, n:(n+numNodesPerHem-1)) = read_niml_binary(getOrigTimeSeriesFilename(subj, hem, experiment));

        n = n + numNodesPerHem;
    end
end
%% END FUNCTION loadOrigTimeseries


%% FUNCTION loadWarpedTimeseries
function T = loadWarpedTimeseries(subj, test_experiment, train_experiment)
    T = zeros(tLen, numNodes);
    
    n = 1;
    for hemNum = 1:nhem
        hem = hems{hemNum};
        T(:, n:(n+numNodesPerHem-1)) = read_niml_binary(getWarpedTimeSeriesFilename(subj, hem, test_experiment, train_experiment));
        
        n = n + numNodesPerHem;
    end
end
%% END FUNCTION loadWarpedTimeseries

end
%% END FUNCTION computeInterSubjectCorrelationResults

