function script_apply_alignment(alignDir, train_experiment, test_experiment)

% Initialize the subjects and hemispheres cell arrays
subjects = {'cb', 'dm', 'hj', 'kd', 'kl', 'mh', 'ph', 'rb', 'se', 'sm'};
hems = {'lh', 'rh'};
outDir = '../surf_data/connnorm_out_data/';
inpSubjDir = '../surf_data/SUBJECTS_DIR/';
saveSubjDir = '../surf_data/SUBJECTS_DIR/';
mmRes = 2;
fwhm = 0;
fs_surf = 'sphere.reg';

train_volDataDir = ['../surf_data/', train_experiment, '/'];
train_suffix = '+orig';
train_extensions = 'BRIK.gz';

test_volDataDir = ['../surf_data/', test_experiment, '/'];
test_suffix = '+orig';
test_extensions = 'BRIK.gz';

preprocess_movie = 0;

% Make the meshes used for aligning test data
make_warped_meshes(subjects, hems, train_experiment, saveSubjDir, outDir, mmRes, fs_surf, alignDir);

% align training data set using training warp
funcnorm_prepare_wrapper(subjects, hems, train_experiment, train_volDataDir, train_suffix, ...
    train_extensions, inpSubjDir, saveSubjDir, preprocess_movie, mmRes, fwhm, ...
    fs_surf, outDir, '', train_experiment, alignDir);

% compute inter-subject correlation results on training set
computeInterSubjectCorrelationResults(subjects, hems, outDir, train_experiment, train_experiment, mmRes, fwhm, fs_surf, alignDir);


% align test data set using training warp
funcnorm_prepare_wrapper(subjects, hems, test_experiment, test_volDataDir, test_suffix, ...
    test_extensions, inpSubjDir, saveSubjDir, preprocess_movie, mmRes, fwhm, ...
    fs_surf, outDir, '', train_experiment, alignDir);

% compute inter-subject correlation results
computeInterSubjectCorrelationResults(subjects, hems, outDir, train_experiment, test_experiment, mmRes, fwhm, fs_surf, alignDir);
