#!/bin/bash

# SLURM job information
# SLURM Job Name- make sure to use the model name (output folder you specify in wrapper!)
#SBATCH -J td_prepro_epicoreg

# Walltime requested
#SBATCH -t 12:00:00

# Number of cores
#SBATCH -n 2

# Memory per core (default = 4GB):
#SBATCH --mem-per-cpu=12G

# Define Oscar partition (or condo) to use
#SBATCH --account=carney-tdesroch-condo

# Provide index values (TASK IDs)
# testing a subject

#SBATCH --array=33,35,36,37,38

# SLURM output and error file names
# Use '%x' for Job Name,'%A' for array-job ID, '%J' for job ID and '%a' for task ID`
#SBATCH -e %x-%A_%a.err
#SBATCH -o %x-%A_%a.out

# Notify user
#SBATCH --mail-user=nadira_yusif_rodriguez@brown.edu
#SBATCH --mail-type=END

# LOAD MODULES AND DEFINE ENVIRONMENT VARIABLES
# get initial time
t0=$(date +%s)
# load modules
module load spm
module load matlab

SESS=`printf "%03d" $SLURM_ARRAY_TASK_ID`

# set environment variables
# paths for input directories
export PREPRO=batchcoregepi
export MONKEY=jesse
export PATH_IN=/gpfs/data/tdesroch/monkey
export PATH_IN_SCRIPTS=/gpfs_home/nyusifro/Github/prepro_monkey
export PATH_IN_BOLD=$PATH_IN/$MONKEY/$SESS/bold
export PATH_IN_ANAT=$PATH_IN/$MONKEY/$SESS/anatomy

# paths for working directories
export PATH_TMP=/tmp/$SESS
export PATH_TMP_SCRIPTS=$PATH_TMP/scripts
export PATH_TMP_SESS=/tmp/$MONKEY/$SESS
export PATH_TMP_ANAT=$PATH_TMP_SESS/anatomy/

# output paths
# use subject's own folder and model name which is SLURM_JOB_NAME
# $SLURM_ARRAY_TASK_ID here is the actual subject number
# $SLURM_JOB_NAME here is the actual name of the model folder
# Removed slurm job name here
export PATH_OUT=/gpfs/data/tdesroch/monkey/$MONKEY/$SESS/bold/

# TRANSFER FILES TO MEMORY (/tmp)
# get time at beginning of transfer step
t1=$(date +%s)
# make temporary directories in shared memory
mkdir $PATH_TMP
mkdir $PATH_TMP_SCRIPTS
mkdir -p $PATH_TMP_SESS
mkdir -p $PATH_TMP_ANAT
# copy input files from gpfs to temporary directory in shared memory
rsync -az $PATH_IN_SCRIPTS/ $PATH_TMP_SCRIPTS
rsync -az $PATH_IN_BOLD/ $PATH_TMP_SESS/bold
rsync -az $PATH_IN_ANAT/ $PATH_TMP_ANAT

# CONDUCT ANALYSIS IN MATLAB
# get time at beginning of analysis step
t2=$(date +%s)

cd $PATH_TMP_SCRIPTS
matlab-threaded -nodisplay -r "td_prepro_coregbatch('$PATH_TMP', '$PATH_TMP_SESS', '$PATH_TMP_ANAT'), exit"

# TRANSFER OUTPUT TO GPFS
# get time at beginning of output transfer step
t3=$(date +%s)
# copy output files from shared memory to gpfs- tried adding a slash bc don't want the whole folder
rsync -az $PATH_TMP_SESS/bold/ $PATH_OUT

# ADDING STEP TO CHANGE FILEPATHS- NOT CURRENTLY RUNNING
# this is now on /gpfs so want to make sure all is setup to do this
# load modules
module load spm
# assuming this script is run from git dir with all other scripts, otherwise need to cd to proper dir
matlab-threaded -nodisplay -r "update_paths_v01('$SESS','$PATH_TMP', '$MONKEY', '$PREPRO'), exit"


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


