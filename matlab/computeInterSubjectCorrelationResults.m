function computeInterSubjectCorrelationResults(subjects, hems, outDir, train_experiment, test_experiment, mmRes, fwhm, fs_surf, alignDir)
% FUNCTION computeInterSubjectCorrelationResults(subjects, hems, outDir, ...
%   train_experiment, test_experiment, mmRes, fwhm, fs_surf [, alignDir])

sep = sprintf('\t');
logFile = '';
incTime = 0;

if ~exist(alignDir) || isempty(alignDir)
	alignDir = [outDir, 'alignment_', train_experiment];
end

if outDir(end) ~= '/'
    outDir = [outDir, '/'];
end

resultsDir = [alignDir, '/results/'];
if isdir(resultsDir) == 0
    mkdir(resultsDir);
end

nsubj = length(subjects);
nhem = length(hems);

tLen = getNumTimepoints(getOrigTimeSeriesFilename(subjects{1}, hems{1}, test_experiment));
numNodesPerHem = getNumNodes(getOrigTimeSeriesFilename(subjects{1}, hems{1}, test_experiment));
numNodes = nhem*numNodesPerHem;

% First, compute the anatomical correlation results
displayLogItem('Computing anatomically-aligned results...', logFile, incTime);
anat_corrs = zeros(nsubj, numNodes);
atl = zeros(tLen, numNodes);
for subjNum = 2:nsubj
    subj = subjects{subjNum};
    
    displayLogItem([sep, 'Loading data for subject #', num2str(subjNum), ' (', subj, ')'], logFile, incTime);
    atl = atl + loadOrigTimeseries(subj, test_experiment);
end

for subjNum = 1:nsubj
    subj = subjects{subjNum};
    
    T = loadOrigTimeseries(subj, test_experiment);
    if subjNum > 1
        atl = atl - T;
    end

    displayLogItem([sep, 'Computing inter-subject correlation results of subject #', num2str(subjNum), ' (', subj, ')'], logFile, incTime);
    for n = 1:numNodes
        magT = sqrt(sum(atl(:, n).^2))*sqrt(sum(T(:, n).^2));
        if magT > 1e-8
            anat_corrs(subjNum, n) = T(:, n)'*atl(:, n)/magT;
        end
    end

    if subjNum < nsubj
        atl = atl + T;
    end
end
% Get the average inter-subject correlation at each node
anat_corrs_avg = sum(anat_corrs)/nsubj;     

% Get the total average inter-subject correlation of each hemisphere
anat_totalAvgCorrPerHem = zeros(1, nhem);
n = 1;
for hemNum = 1:nhem
    hem = hems{hemNum};
    anat_totalAvgCorrPerHem(hemNum) = sum(anat_corrs_avg(n:(n+numNodesPerHem-1)))/numNodesPerHem;
    
    outStr = [resultsDir, test_experiment, '_anatAlign_avg_intersubjectcorrelation_', hem, '.dset'];
    fp = fopen(outStr, 'w');
    fprintf(fp, '%f\n', anat_corrs_avg(n:(n+numNodesPerHem-1)));
    fclose(fp);
    
    displayLogItem(['Anatomically-aligned inter-subject correlation results (', hem, ' hemisphere) saved to:  ', outStr]);
    
    n = n + numNodesPerHem;
end
% Get the total average inter-subject correlation across all hemispheres
anat_totalAvgCorr = sum(anat_totalAvgCorrPerHem)*numNodesPerHem/numNodes;


% Now, compute the updated warped results
displayLogItem('Computing functionally-aligned results...', logFile, incTime);
func_corrs = zeros(nsubj, numNodes);
atl = zeros(tLen, numNodes);
for subjNum = 2:nsubj
    subj = subjects{subjNum};
    
    displayLogItem([sep, 'Loading data for subject #', num2str(subjNum), ' (', subj, ')'], logFile, incTime);
    atl = atl + loadWarpedTimeseries(subj, test_experiment, train_experiment);
end

for subjNum = 1:nsubj
    subj = subjects{subjNum};
    
    displayLogItem([sep, 'Computing inter-subject correlation results of subject #', num2str(subjNum), ' (', subj, ')'], logFile, incTime);
    T = loadWarpedTimeseries(subj, test_experiment, train_experiment);
    if subjNum > 1
        atl = atl - T;
    end
    
    for n = 1:numNodes
        magT = sqrt(sum(atl(:, n).^2))*sqrt(sum(T(:, n).^2));
        if magT > 1e-8
            func_corrs(subjNum, n) = T(:, n)'*atl(:, n)/magT;
        end
    end
    
    if subjNum < nsubj
        atl = atl + T;
    end
end
% Get the average inter-subject correlation at each node
func_corrs_avg = sum(func_corrs)/nsubj;     

% Get the total average inter-subject correlation of each hemisphere
func_totalAvgCorrPerHem = zeros(1, nhem);
n = 1;
for hemNum = 1:nhem
    hem = hems{hemNum};
    func_totalAvgCorrPerHem(hemNum) = sum(func_corrs_avg(n:(n+numNodesPerHem-1)))/numNodesPerHem;
    
    outStr = [resultsDir, test_experiment, '_funcAlign_warpedby_', train_experiment, '_avg_intersubjectcorrelation_', hem, '.dset'];
    fp = fopen(outStr, 'w');
    fprintf(fp, '%f\n', func_corrs_avg(n:(n+numNodesPerHem-1)));
    fclose(fp);
    
    displayLogItem(['Functionally-aligned inter-subject correlation results (', hem, ' hemisphere) saved to:  ', outStr]);
    
    n = n + numNodesPerHem;
end
% Get the total average inter-subject correlation across all hemispheres
func_totalAvgCorr = sum(func_totalAvgCorrPerHem)*numNodesPerHem/numNodes;


displayLogItem(sprintf('\n'), logFile, incTime);
displayLogItem(['Anatomical Alignment Results (on ', test_experiment, '):'], logFile, incTime);
    displayLogItem([sep, 'Average inter-subject correlation (per hemisphere):  '], logFile, incTime);
        for hemNum = 1:nhem
            hem = hems{hemNum};
            displayLogItem([sep, sep, hem, ':  ', num2str(anat_totalAvgCorrPerHem(hemNum))], logFile, incTime);
        end
    displayLogItem([sep, 'Total average inter-subject correlation:  ', num2str(anat_totalAvgCorr)], logFile, incTime);

displayLogItem(sprintf('\n'), logFile, incTime);
displayLogItem(['Functional Alignment Results (on ', test_experiment, '):'], logFile, incTime);
    displayLogItem([sep, 'Average inter-subject correlation (per hemisphere):  '], logFile, incTime);
        for hemNum = 1:nhem
            hem = hems{hemNum};
            displayLogItem([sep, sep, hem, ':  ', num2str(func_totalAvgCorrPerHem(hemNum))], logFile, incTime);
        end
    displayLogItem([sep, 'Total average inter-subject correlation:  ', num2str(func_totalAvgCorr)], logFile, incTime);





    
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