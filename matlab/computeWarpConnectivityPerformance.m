function computeWarpConnectivityPerformance(atlasDir, testExp, passNum, zeroCorrected)
% FUNCTION computeWarpConnectivityPerformance(atlasDir, testExp, passNum, zeroCorrected)

isAnat = 0;
if strcmp(atlasDir, 'anat/')
    isAnat = 1;
    passNum = 0;
    zeroCorrected = '';

    fileSuffix = '';
else
    if zeroCorrected == 1
        corrected = '_corrected';
    else
        corrected = '';
    end

    fileSuffix = ['_pass_', num2str(passNum), corrected];
end

if strcmp(testExp, 'raidersP2')
    tLen = 1112;
    testExpStr = 'raiders/r2';
    tFileName = 'noglobal_on_std_sphere.reg.float_x86_64_bin';
elseif strcmp(testExp, 'raidersP1')
    tLen = 1101;
    testExpStr = 'raiders/r1';
    tFileName = 'noglobal_on_std_sphere.reg.float_x86_64_bin';
elseif strcmp(testExp, 'monkeydog')
    tLen = 1536;
    testExpStr = 'monkeydog/r1';
    tFileName = 'std_dt_F1_on_std_sphere.reg.float_x86_64_bin';
else
    error('unknown test exp');
end


subjects = {'cb', 'dm', 'hj', 'kd', 'kl', 'mh', 'ph', 'rb', 'se', 'sm'};
nsubj = length(subjects);
hems = {'lh', 'rh'};
nhem = length(hems);

numNodesPerHem = 36002;
numNodes = nhem*numNodesPerHem;

loopStep = 500;

outDir = [atlasDir, 'results_connectivity_correlation/'];
if isAnat == 0
    outDir = [outDir, 'pass_', num2str(passNum), corrected, '/'];
end
if isdir(outDir) == 0
    mkdir(outDir);
end

avgCorrs = zeros(nhem, nsubj);

% First check if this has already been computed
completedAlready = 0;
doneCount = 0;
for j = 1:nsubj
    subj = subjects{j};
    
    for k = 1:nhem
        hem = hems{k};
        
        e = exist([outDir, 'connectivitycorr5_', hem, '_', testExp, '_subj_', subj, fileSuffix, '.bin']);
        if e == 2
            % Then it exists
            fp = fopen([outDir, 'connectivitycorr5_', hem, '_', testExp, '_subj_', subj, fileSuffix, '.bin'], 'rb');
            avgCorrs(k, j) = sum(fread(fp, [numNodesPerHem, 1], 'float32'))/numNodesPerHem;
            fclose(fp);
            doneCount = doneCount + 1;
        end
    end
end

if doneCount == nhem*nsubj
    completedAlready = 1;
end

if completedAlready == 1
    avgCorrs
    sum(avgCorrs')/nsubj
    return;
end

fp = fopen('../fmridata/standards/std_2mm_36002nodes/cartesianCoords.bin', 'rb');
cartCoords = fread(fp, [3, 36002], 'float32');fclose(fp);
for j = 1:36002; cartCoords(:, j) = cartCoords(:, j)/sqrt(sum(cartCoords(:, j).^2));end;
fp = fopen('../fmridata/standards/std_2mm_36002nodes/nearest_1mm_144002node.bin', 'rb');
nearestNodes = fread(fp, [1, 36002], 'int32');fclose(fp);

fp = fopen('../fmridata/standards/std_1mm_144002nodes/cartesianCoords.bin', 'rb');
cartCoords_hires = fread(fp, [3, 144002], 'float32');fclose(fp);
for j = 1:144002; cartCoords_hires(:, j) = cartCoords_hires(:, j)/sqrt(sum(cartCoords_hires(:, j).^2));end;
fp = fopen('../fmridata/standards/std_1mm_144002nodes/neighbors.bin', 'rb');
neighbors_hires = fread(fp, [6, 144002], 'int32');fclose(fp);

disp('Computing functional atlas...');pause(0.05);
for subjNum = 1:nsubj
    subj = subjects{subjNum};
    
    if isAnat == 0
        fp = fopen([atlasDir, 'warp_subj_', subj, '_to_atlas_pass_', num2str(passNum), corrected], 'rb');
        warp{1} = fread(fp, [3, numNodesPerHem], 'float32');
        warp{2} = fread(fp, [3, numNodesPerHem], 'float32');
        fclose(fp);
    end
    
    T = zeros(tLen, numNodes);
    for hemNum = 1:nhem
        hem = hems{hemNum};
        
        kmin = numNodesPerHem*(hemNum - 1) + 1;
        kmax = kmin + numNodesPerHem - 1;

        if isAnat == 0
            fp = fopen(['/data/subjects_out/', subj, '/funcnorm/', testExpStr, '/', hem, '/std_sphere.reg_144002/sm0/', tFileName], 'rb');
            T_hires = fread(fp, [tLen, 144002], 'float32');fclose(fp);

            T(:, kmin:kmax) = computeInterpOnSphere_hires(T_hires, cartCoords, warp{hemNum}, nearestNodes, cartCoords_hires, neighbors_hires);
        else
            fp = fopen(['../fmridata/', testExp, '_', subj, '_reg_', hem, '/timeseries.bin'], 'rb');
            T(:, kmin:kmax) = fread(fp, [tLen, numNodesPerHem], 'float32');fclose(fp);
        end        
    end

    % Normalize time-series
    T = normalizeTimeseries(T, tLen, numNodes);
    
    % Write out the normalized timeseries
    fp = fopen([outDir, 'subjwarped_', subj, '_', testExp, fileSuffix, '.bin'], 'wb');
    fwrite(fp, T, 'float32');fclose(fp);
        
    disp(['Completed subject #', num2str(subjNum), ' out of ', num2str(nsubj)]);pause(0.05);
end
disp('Completed computing functional atlas!');pause(0.05);

% Clear all connectivity corr files
for subjNum = 1:nsubj
    subj = subjects{subjNum};    
    for hemNum = 1:nhem
        hem = hems{hemNum};
        fp = fopen([outDir, 'connectivitycorr5_', hem, '_', testExp, '_subj_', subj, fileSuffix, '.bin'], 'wb');fclose(fp);
    end
end

corrs = zeros(nsubj, numNodes);
numLoops = ceil(numNodes/loopStep);
for k = 1:loopStep:numNodes
    disp(['Working on loop #', num2str((k-1)/loopStep + 1), ' out of ', num2str(numLoops), '...']);pause(0.05);
    
    kMax = min(k+loopStep-1, numNodes);
    kInd = k:kMax;    
    currNumNodes = length(kInd);

    atl = zeros(numNodes, currNumNodes);    
    subj_means = zeros(nsubj, currNumNodes);
    subj_sds = zeros(nsubj, currNumNodes);

    for subjNum = 1:nsubj
        c{subjNum} = zeros(numNodes, currNumNodes);
    end

    disp(['Generating leave-one-out-atlas...']);pause(0.05);
    for subjNum = 2:nsubj
        subj = subjects{subjNum};
        
        fp = fopen([outDir, 'subjwarped_', subj, '_', testExp, fileSuffix, '.bin'], 'rb');
        T = fread(fp, [tLen, numNodes], 'float32');fclose(fp);
        
        c{subjNum} = T'*T(:, kInd);
        subj_means(subjNum, :) = sum(c{subjNum})/numNodes;
        subj_sds(subjNum, :) = sqrt(sum(c{subjNum}.^2) - numNodes*subj_means(subjNum, :).^2);
        
        atl = atl + c{subjNum}/(nsubj-1);
    end
    disp(['Completed generating leave-one-out-atlas!']);pause(0.05);
    

    disp(['Computing correlations...']);pause(0.05);
%    atl_means = zeros(1, numNodes);
%    atl_sds = zeros(1, numNodes);
    for subjNum = 1:nsubj
        if subjNum == 1
            % Then just compute correlation
            fp = fopen([outDir, 'subjwarped_', subj, '_', testExp, fileSuffix, '.bin'], 'rb');
            T = fread(fp, [tLen, numNodes], 'float32');fclose(fp);
            
            c{subjNum} = T'*T(:, kInd);
%            subj_means(subjNum, :) = sum(c{subjNum})/numNodes;
%            subj_sds(subjNum, :) = sqrt(sum(c{subjNum}.^2) - numNodes*subj_means(subjNum, :).^2);
        else
            % Remove the subject from the atlas, and add in previous
            % subject
            atl = atl + (c{subjNum-1} - c{subjNum})/(nsubj - 1);
        end
        
        % Compute mean and s.d. of current atlas
%        atl_means = sum(atl)/numNodes;
%        atl_sds = sqrt(sum(atl.^2) - numNodes*atl_means.^2);
        
        sC = c{subjNum};
        
        for n = 1:currNumNodes
%            if atl_sds(n) < 1e-8 || subj_sds(subjNum, n) < 1e-8
%                continue;
%            end

            % corrcoef returns 2 x 2 matrix R
            R = corrcoef(sC(:, n), atl(:, n));
            corrs(subjNum, kInd(n)) = R(2);
%            corrs(subjNum, kInd(n)) = (sC(:, n)'*atl(:, n) - numNodes*atl_means(n)*subj_means(subjNum, n))/(atl_sds(n)*subj_sds(subjNum, n));
        end
    end
    disp(['Completed computing correlations!']);pause(0.05);
    
    disp(['Writing out correlations...']);pause(0.05);
    for subjNum = 1:nsubj
        subj = subjects{subjNum};
        
        kMax2 = k;
        if k < numNodesPerHem
            kMax2 = min(kMax, numNodesPerHem);
            fp = fopen([outDir, 'connectivitycorr5_lh_', testExp, '_subj_', subj, fileSuffix, '.bin'], 'ab');
            fwrite(fp, corrs(subjNum, k:kMax2), 'float32');fclose(fp);
            kMax2 = kMax2 + 1;
        end
        
        if kMax2 < kMax
            fp = fopen([outDir, 'connectivitycorr5_rh_', testExp, '_subj_', subj, fileSuffix, '.bin'], 'ab');
            fwrite(fp, corrs(subjNum, kMax2:kMax), 'float32');fclose(fp);
        end
    end    
    disp(['Completed writing out correlations!']);
    
    disp(['Completed loop #', num2str((k-1)/loopStep + 1), ' out of ', num2str(numLoops), '!']);pause(0.05);
end
avgCorrs = zeros(nhem, nsubj);

for subjNum = 1:nsubj
    subj = subjects{subjNum};
    
    for hemNum = 1:nhem
        hStart = (hemNum-1)*numNodesPerHem + 1;
        hEnd = hStart + numNodesPerHem - 1;
        
        avgCorrs(hemNum, subjNum) = sum(corrs(subjNum, hStart:hEnd)')/numNodesPerHem;
    end
end
avgCorrs

function T = normalizeTimeseries(T, tLen, numNodes)

for j = 1:numNodes
    T(:, j) = T(:, j) - sum(T(:, j))/tLen;
    mag1 = sqrt(sum(T(:, j).^2));
    if mag1 < 1e-8
        mag1 = 0;
    else
        mag1 = 1/mag1;
    end
    
    T(:, j) = mag1*T(:, j);
end