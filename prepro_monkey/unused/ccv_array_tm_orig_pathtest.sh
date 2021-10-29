#!/bin/bash

# SLURM job information
# SLURM Job Name- make sure to use the model name (output folder you specify in wrapper!)
#SBATCH -J glm_para1234_v15_2

# Walltime requested
#SBATCH -t 12:00:00

# Number of cores
#SBATCH -n 2

# Memory per core (default = 4GB):
#SBATCH --mem-per-cpu=12G

# Define Oscar partition (or condo) to use
#SBATCH -A bibs-tdesroch-condo

# Provide index values (TASK IDs)
#SBATCH --array=107
# testing 3 subjects

# SLURM output and error file names
# Use '%x' for Job Name,'%A' for array-job ID, '%J' for job ID and '%a' for task ID`
#SBATCH -e %x-%A_%a.err
#SBATCH -o %x-%A_%a.out

# Notify user
#SBATCH --mail-user=theresa_mckim@brown.edu
#SBATCH --mail-type=END

# LOAD MODULES AND DEFINE ENVIRONMENT VARIABLES
# get initial time
t0=$(date +%s)
# load modules
module load spm

# set environment variables
# paths for input directories
export PATH_IN=/gpfs/data/tdesroch/rmSEQ
export PATH_IN_SCRIPTS=/gpfs_home/tmckim/rmseq_analysis
export PATH_IN_ANAT=$PATH_IN/subjects/$SLURM_ARRAY_TASK_ID/anatomy
export PATH_IN_BEHAV=$PATH_IN/subjects/$SLURM_ARRAY_TASK_ID/behavior
export PATH_IN_BOLD=$PATH_IN/subjects/$SLURM_ARRAY_TASK_ID/bold

# paths for working directories
export PATH_TMP=/tmp/$SLURM_JOB_ID
export PATH_TMP_SCRIPTS=$PATH_TMP/scripts
export PATH_TMP_SUBJECTS=$PATH_TMP/subjects

# output paths
# use subject's own folder and model name which is SLURM_JOB_NAME
# $SLURM_ARRAY_TASK_ID here is the actual subject number
# $SLURM_JOB_NAME here is the actual name of the model folder
# Trying to remove the last part of this path $SLURM_JOB_NAME to see if only one model folder name gets made in the subj level folder, or if it just adds the files without a folder
export PATH_OUT=/gpfs/data/tdesroch/rmSEQ/subjects/$SLURM_ARRAY_TASK_ID/
# Trying to remove the last part of this path /$SLURM_JOB_NAME/$SLURM_ARRAY_TASK_ID to see if folders don't get made in the group level folder
export PATH_OUT_GROUP=/gpfs/data/tdesroch/rmSEQ/GroupAnalysis

# TRANSFER FILES TO MEMORY (/tmp)
# get time at beginning of transfer step
t1=$(date +%s)
# make temporary directories in shared memory
mkdir $PATH_TMP
mkdir $PATH_TMP_SCRIPTS
mkdir -p $PATH_TMP_SUBJECTS/$SLURM_ARRAY_TASK_ID
# copy input files from gpfs to temporary directory in shared memory
rsync -az $PATH_IN_SCRIPTS/ $PATH_TMP_SCRIPTS
rsync -az $PATH_IN_ANAT $PATH_IN_BEHAV $PATH_IN_BOLD $PATH_TMP_SUBJECTS/$SLURM_ARRAY_TASK_ID

# CONDUCT ANALYSIS IN MATLAB
# get time at beginning of analysis step
t2=$(date +%s)

cd $PATH_TMP_SCRIPTS
matlab-threaded -nodisplay -r "rms_GLMwrapper_para1234_v16([true false false], $SLURM_ARRAY_TASK_ID,'$PATH_TMP'), exit"

# TRANSFER OUTPUT TO GPFS
# get time at beginning of output transfer step
t3=$(date +%s)
# create directory structure for output files
# removing the mkdir for pathout
# mkdir $PATH_OUT
mkdir -p $PATH_OUT_GROUP
# copy output files from shared memory to gpfs
rsync -az $PATH_TMP_SUBJECTS/$SLURM_ARRAY_TASK_ID/$SLURM_JOB_NAME $PATH_OUT
rsync -az $PATH_TMP/GroupAnalysis/ $PATH_OUT_GROUP


# ADDING STEP TO CHANGE FILEPATHS
# this is now on /gpfs so want to make sure all is setup to do this
# load modules
module load spm
# assuming this script is run from git dir with all other scripts, otherwise need to cd to proper dir
matlab-threaded -nodisplay -r "rms_updatepaths_v01($SLURM_ARRAY_TASK_ID,'$SLURM_JOB_NAME','$PATH_TMP'), exit"


# REMOVE TEMPORARY FILES FROM MEMORY (/tmp)
# get time at beginning of cleanup
t4=$(date +%s)
# remove files from memory
rm -rf $PATH_TMP

# OUTPUT TIMING INFORMATION (optional, for diagnostic purposes)
dt1=$(echo "$t1 - $t0" | bc)
dt2=$(echo "$t2 - $t1" | bc)
dt3=$(echo "$t3 - $t2" | bc)
dt4=$(echo "$t4 - $t3" | bc)
dt=$(echo "$t4 - $t0" | bc)

printf "Initialization (s): %08.1f\n" $dt1
printf "Transfer input from gpfs (s): %08.1f\n" $dt2
printf "Analysis step (s): %08.1f\n" $dt3
printf "Transfer output to gpfs (s): %08.1f\n" $dt4
printf "Total (s): %08.1f\n" $dt


