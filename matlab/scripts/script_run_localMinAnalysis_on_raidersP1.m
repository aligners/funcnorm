function script_run_localMinAnalysis_on_raidersP1(lambda, maxResolution)
% FUNCTION script_run_withinSubjReg_on_raidersP1
% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

if nargin == 1
    maxResolution = 3;
end
% Initialize the subjects and hemispheres cell arrays
subjects = {'cb', 'dm', 'hj', 'kd', 'kl', 'mh', 'ph', 'rb', 'se', 'sm'};
hems = {'lh', 'rh'};

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

mmRes = 2;
fwhm = 0;
fs_surf = 'sphere.reg';
outDir = '../surf_data/connnorm_out_data/';
lambda_areal = lambda;
lambda_metric = lambda;
logFile = 'connnorm_localMinAnalysis_log.txt';


script_connnorm_localMinAnalysis(subjects, hems, train_experiment, mmRes, fwhm, fs_surf, outDir, ...
   lambda_areal, lambda_metric, maxResolution, logFile);
