function script_run_connnorm_on_raidersP1(lambda, tLen)
% FUNCTION script_run_connnorm_on_raidersP1(lambda, tLen)
% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

% Initialize the subjects and hemispheres cell arrays
subjects = {'cb', 'dm', 'hj', 'kd', 'kl', 'mh', 'ph', 'rb', 'se', 'sm'};
hems = {'lh', 'rh'};
pctEnergy = 95;

if nargin == 2
	hasTLen = 1;
else
	hasTLen = 0;
end

% Set the training data (dataset used to define the functional correspondence)
train_volDataDir = '../surf_data/raidersP1/';
train_experiment = 'raidersP1';
train_suffix = '+orig';
train_extensions = 'BRIK.gz';

% Set the test data (independent dataset that will be used for evaluating the alignment
% derived on the training set)
test_volDataDir = '../surf_data/raidersP2/';
test_experiment = 'raidersP2';
test_suffix = '+orig';
test_extensions = 'BRIK.gz';

inpSubjDir = '../surf_data/SUBJECTS_DIR/';
saveSubjDir = '../surf_data/SUBJECTS_DIR/';

preprocess_movie = 0;
mmRes = 2;
fwhm = 0;
fs_surf = 'sphere.reg';
outDir = '../surf_data/connnorm_out_data/';
lambda_areal = lambda;
lambda_metric = lambda;
maxResolution = 3;
numPasses = 3;
logFile = 'funcnorm_log.txt';

if hasTLen == 1
	connnorm(subjects, hems, train_experiment, mmRes, fwhm, fs_surf, outDir, ...
		lambda_areal, lambda_metric, pctEnergy, maxResolution, numPasses, logFile, tLen);	
else
	connnorm(subjects, hems, train_experiment, mmRes, fwhm, fs_surf, outDir, ...
		lambda_areal, lambda_metric, pctEnergy, maxResolution, numPasses, logFile);
end

% Make the meshes used for aligning test data
%make_warped_meshes(subjects, hems, train_experiment, saveSubjDir, outDir, mmRes, fs_surf);

% align training data set using training warp
%funcnorm_prepare_wrapper(subjects, hems, train_experiment, train_volDataDir, train_suffix, ...
%    train_extensions, inpSubjDir, saveSubjDir, preprocess_movie, mmRes, fwhm, ...
%    fs_surf, outDir, '', train_experiment);

% compute inter-subject correlation results on training set
%computeInterSubjectCorrelationResults(subjects, hems, outDir, train_experiment, train_experiment, mmRes, fwhm, fs_surf);


% align test data set using training warp
%funcnorm_prepare_wrapper(subjects, hems, test_experiment, test_volDataDir, test_suffix, ...
%    test_extensions, inpSubjDir, saveSubjDir, preprocess_movie, mmRes, fwhm, ...
%    fs_surf, outDir, '', train_experiment);

% compute inter-subject correlation results
%computeInterSubjectCorrelationResults(subjects, hems, outDir, train_experiment, test_experiment, mmRes, fwhm, fs_surf);

