function demo_funcnorm
%function demo_funcnorm
% Runs all steps of the funcnorm demo

% Copyright (c) 2008 Bryan Conroy, Jim Haxby, Peter Ramadge, Mert Rory Sabuncu, and Ben Singer

% Initialize the subjects and hemispheres cell arrays
subjects = {'cb', 'dm', 'hj', 'kd', 'mh', 'rb'};
hems = {'lh', 'rh'};

% Set the training data (dataset used to define the functional correspondence)
train_volDataDir = '../sample_data/movie/';
train_experiment = 'movie5min_P1';
train_suffix = '+orig';
train_extensions = 'BRIK.gz';

% Set the test data (independent dataset that will be used for evaluating the alignment
% derived on the training set)
test_volDataDir = '../sample_data/movie/';
test_experiment = 'movie5min_P2';
test_suffix = '+orig';
test_extensions = 'BRIK.gz';

inpSubjDir = '../sample_data/SUBJECTS_DIR/';
saveSubjDir = '../sample_results/subjects_dir/';
preprocess_movie = 1;
mmRes = 2;
fwhm = 0;
fs_surf = 'sphere.reg';
outDir = '../sample_results/outputdir/';
lambda_areal = 30;
lambda_metric = 30;
maxResolution = 3;
numPasses = 3;
logFile = 'funcnorm_log.txt';

inpUnixPath=getenv('PATH');
if ismac
    % On the Mac your ~/.bashrc, ~/.cshrc, etc are ignored when Matlab
    % is launched by double-clicking on its icon, so setting up your PATH
    % in login files won't help. Prepend the paths to your particular versions of
    % FreeSurfer and Afni here, unless you've figured out the official method.
    % 
    % Here's mine, you'll want to change it to yours if you are on a Mac:
    inpUnixPath=['/Applications/freesurfer/bin:/Users/bdsinger/abin:',inpUnixPath];
end

% Run the alignment on the training data -- the dataset will first be prepared by funcnorm_prepare.sh
funcnorm(subjects, hems, train_experiment, train_volDataDir, train_suffix, ...
    train_extensions, inpSubjDir, saveSubjDir, preprocess_movie, mmRes, fwhm, fs_surf, ...
    outDir, lambda_areal, lambda_metric, maxResolution, numPasses, logFile, ...
    '', inpUnixPath);

% Prepare the test dataset (movie P2)
funcnorm_prepare_wrapper(subjects, hems, test_experiment, test_volDataDir, ...
    test_suffix, test_extensions, inpSubjDir, saveSubjDir, preprocess_movie, mmRes, ...
    fwhm, fs_surf, outDir, '', inpUnixPath);

% Make the meshes used for aligning test data
make_warped_meshes(subjects, hems, train_experiment, saveSubjDir, outDir, mmRes, fs_surf);

% align test data set using training warp
funcnorm_prepare_wrapper(subjects, hems, test_experiment, test_volDataDir, test_suffix, ...
    test_extensions, inpSubjDir, saveSubjDir, preprocess_movie, mmRes, fwhm, ...
    fs_surf, outDir, '', train_experiment, inpUnixPath);

% compute inter-subject correlation results
computeInterSubjectCorrelationResults(subjects, hems, outDir, train_experiment, test_experiment, mmRes, fwhm, fs_surf);

