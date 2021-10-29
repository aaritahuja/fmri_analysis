function task_first_levels(MRI_directory, subjids, task_runs, varargin)
%MRI_directory is the directory where all the MRI data is stored
%('/shared/lab/projects/analysis/aarit/MRI_data')
%subjids is a matrix of all the subjects to be looped over
%good_runs is a cell array of all of the task runs for each subject
%example: good_runs = {[3:7, 9:15], 3:14}
%localizer_runs is a cell array of all of the task runs for each subject
%example: localizer_runs = {1:2, 1:2, 1, 1:2}

%For the first batch of 12 human subjects the variables are
%subjids = [2:7, 10:14, 16];
%task_runs = {[3:7, 9:15], 3:14, 3:13, 3:14, 3:10, 3:14, 2:13, 2:13, 2:13, 2:13, [2:7, 9:10], 2:13}
%localizer_runs = {1:2, 1:2, 1:2, 1:2, 1:2, 1:2, 1, 1, 1, 1, 1, 1};

%put in something that checks that the length of subjruns and subjids is
%the same

for i = 1:length(subjids)
    subjid = subjids(i);
    task_run = cell2mat(task_runs(i));
    %localizer_run = cell2mat(localizer_runs(i));
    SPM_path = task_design_matrix(MRI_directory, subjid, task_run);
    contrast_generator(SPM_path, 'simulation', '')
    contrast_generator(SPM_path, 'perception', '')
    contrast_generator(SPM_path, 'counting', '')
    contrast_generator(SPM_path, 'simulation', 'counting')
    contrast_generator(SPM_path, 'perception', 'counting')
    contrast_generator(SPM_path, 'mvt', '')
    %     contrast_generator(SPM_path, 'simulation*bf', 'counting', 'simulationxuncertainty')
end
