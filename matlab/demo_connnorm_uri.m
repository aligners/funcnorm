function demo_connnorm_uri(trainExpt, testExpt)

% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

% Initialize the subjects and hemispheres cell arrays
subjects = {'bc_uri', 'ch_uri', 'ls_uri'};
hems = {'lh', 'rh'};
if strcmp(trainExpt, 'darkness')
    pctEnergy = 95;
elseif strcmp(trainExpt, 'story2')
    pctEnergy = 95;
else
    error(['unknown expt ', trainExpt]);
end

% Set the training data (dataset used to define the functional correspondence)
if strcmp(trainExpt, 'darkness')
	train_volDataDir = '/home/bconroy/uri_data/darkness/';
	train_experiment = 'darkness';
elseif strcmp(trainExpt, 'story2')
    train_volDataDir = '/home/bconroy/uri_data/story2/';
    train_experiment = 'story2';
else
    error(['unknown expt ', trainExpt]);
end
train_suffix = '+orig';
train_extensions = 'BRIK.gz';

% Set the test data (independent dataset that will be used for evaluating the alignment
% derived on the training set)
if strcmp(testExpt, 'darkness')
	test_volDataDir = '/home/bconroy/uri_data/darkness/';
	test_experiment = 'darkness';
elseif strcmp(testExpt, 'story2')
    test_volDataDir = '/home/bconroy/uri_data/story2/';
    test_experiment = 'story2';
else
	error(['unknown expt:  ', testExpt]);
end
test_suffix = '+orig';
test_extensions = 'BRIK.gz';

inpSubjDir = '/data/subjects_in/';
saveSubjDir = '/home/bconroy/uri_data/sample_results/subjects_dir/';
preprocess_movie = 0;
mmRes = 2;
fwhm = 0;
fs_surf = 'sphere.reg';
outDir = '/home/bconroy/uri_data/sample_results/';
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
connnorm(subjects, hems, train_experiment, mmRes, fwhm, fs_surf, outDir, ...
   lambda_areal, lambda_metric, pctEnergy, maxResolution, numPasses, logFile);

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

