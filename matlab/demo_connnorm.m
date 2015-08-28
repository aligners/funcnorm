function demo_connnorm(expt)

% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

% Initialize the subjects and hemispheres cell arrays
if strcmp(expt, 'movie5min_P1')
	subjects = {'cb', 'dm', 'hj', 'kd', 'mh', 'rb'};
elseif strcmp(expt, 'raidersP1')
	subjects = {'cb', 'dm', 'hj', 'kd', 'kl', 'mh', 'ph', 'rb', 'se', 'sm'};
else
	error(['unknown expt:  ', expt]);
end

hems = {'lh', 'rh'};
%pctEnergy = 95;
pctEnergy = 80;

% Set the training data (dataset used to define the functional correspondence)
if strcmp(expt, 'movie5min_P1')
	train_volDataDir = '../sample_data/movie/';
	train_experiment = 'movie5min_P1';
elseif strcmp(expt, 'raidersP1')
	train_volDataDir = '../sample_data/raidersP1/';
	train_experiment = 'raidersP1';
else
	error(['unknown expt:  ', expt]);
end
train_suffix = '+orig';
train_extensions = 'BRIK.gz';

% Set the test data (independent dataset that will be used for evaluating the alignment
% derived on the training set)
if strcmp(expt, 'movie5min_P1')
	test_volDataDir = '../sample_data/movie/';
	test_experiment = 'movie5min_P2';
elseif strcmp(expt, 'raidersP1')
	test_volDataDir = '../sample_data/raidersP2/';
	test_experiment = 'raidersP2';
else
	error(['unknown expt:  ', expt]);
end
test_suffix = '+orig';
test_extensions = 'BRIK.gz';

inpSubjDir = '../sample_data/SUBJECTS_DIR/';
saveSubjDir = '../sample_results/subjects_dir/';
preprocess_movie = 1;
mmRes = 2;
fwhm = 0;
fs_surf = 'sphere.reg';
outDir = '../sample_results/';
lambda_areal = 4.17;
lambda_metric = 4.17;
maxResolution = 3;
numPasses = 3;
logFile = 'connorm_log.txt';

% Run the alignment on the training data -- the dataset will first be prepared by funcnorm_prepare.sh
%connnorm(subjects, hems, train_experiment, train_volDataDir, train_suffix, ...
%    train_extensions, inpSubjDir, saveSubjDir, preprocess_movie, mmRes, fwhm, fs_surf, ...
%    outDir, lambda_areal, lambda_metric, pctEnergy, maxResolution, numPasses, logFile);

% If we didn't need to run funcnorm_prepare, then use the line below
connnorm('subjects',subjects, 'hems',hems, 'train_experiment',train_experiment, 'mmRes',mmRes, ...
            'fwhm',fwhm, 'fs_surf',fs_surf, 'outDir',outDir, 'lambda_areal',lambda_areal, ...
            'lambda_metric',lambda_metric, 'pctEnergy',pctEnergy, 'maxResolution',maxResolution, ...
            'numPasses',numPasses, 'logFile',logFile);

% Prepare the test dataset (movie P2)
funcnorm_prepare_wrapper(subjects, hems, test_experiment, test_volDataDir, ...
    test_suffix, test_extensions, inpSubjDir, saveSubjDir, preprocess_movie, mmRes, ...
    fwhm, fs_surf, outDir);

% Make the meshes used for aligning test data
make_warped_meshes(subjects, hems, train_experiment, saveSubjDir, outDir, mmRes, fs_surf);

% align test data set using training warp
funcnorm_prepare_wrapper(subjects, hems, test_experiment, test_volDataDir, test_suffix, ...
    test_extensions, inpSubjDir, saveSubjDir, preprocess_movie, mmRes, fwhm, ...
    fs_surf, outDir, '', train_experiment);

% compute inter-subject correlation results
computeInterSubjectCorrelationResults(subjects, hems, outDir, train_experiment, test_experiment, mmRes, fwhm, fs_surf);

