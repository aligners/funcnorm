#!/bin/bash

# This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
# Please see AUTHORS and LICENSE file in the project root directory

# where sample data lives
sample_data_dir="../sample_data"
sample_results_dir="../sample_results"

# save original SUBJECTS_DIR and restore at end of script
orig_subj_dir=$SUBJECTS_DIR
export SUBJECTS_DIR=".$sample_data_dir/SUBJECTS_DIR"

# local constants for this preparation
inputdir="$sample_data_dir/movie"
experiment="movie5min_P1"
suffix="+orig"
subjects_dir="$sample_results_dir/subjects_dir"
outputdir="$sample_results_dir/outputdir"
fs_surf="sphere.reg"
mm=2
fwhm=0
verbosity=2

for sub in cb dm hj kd mh rb; do
	for hem in lh rh; do
		for part in P1 P2; do
			echo "Preparing movie $part data for subject $sub, $hem hemisphere"
			./funcnorm_prepare.sh \
				-sub $sub \
				-hem $hem \
				-inputdir $inputdir
				-experiment $experiment \
				-suffix $suffix \
				-subjects_dir $subjects_dir \
				-outputdir $outputdir \
				-fs_surf $fs_surf \
				-mm $mm \
				-fwhm $fwm  \
				-preprocess_movie \
				-verbosity $verbosity
				
				# check that it ran okay
				if [ $? -ne 0 ]; then
					echo "The script funcnorm_prepare.sh failed. " >&2
					export SUBJECTS_DIR=$orig_subj_dir
					exit -1
				fi
		done
	done
done

# restore SUBJECTS_DIR
export SUBJECTS_DIR=$orig_subj_dir

echo "Data has been prepared. Now start Matlab and run demo_funcnorm.m"
