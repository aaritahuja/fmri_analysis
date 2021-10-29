function localizer_first_levels(MRI_directory, subjids, localizer_runs, sessions, hrf_type, preprocessing_package)
%MRI_directory is the directory where all the MRI data is stored
%('/shared/lab/projects/analysis/aarit/MRI_data')
%subjids is a matrix of all the subjects to be looped over
%localizer_runs is a cell array of all of the localizer runs for each subject
%example: localizer_runs = {1:2, 1:2, 1, 1:2}

%put in something that checks that the length of subjruns and subjids is
%the same
if isnumeric(subjids)
    iterate = subjids;
else
    iterate = sessions;
end

for i = 1:length(iterate)
    if isnumeric(subjids)
        subjid = subjids(i);
    else
        subjid = subjids;
    end
    l_runs = localizer_runs{i};
    session = strcat('sub-glenn', sprintf('%03d', sessions(i)), '/ses-01');
    SPM_path = localizer_design_matrix(MRI_directory, subjid, l_runs, session, hrf_type, preprocessing_package);
    contrast_generator(SPM_path, 'motion', 'static')
    %contrast_generator(SPM_path, 'static', 'motion')
end

%contrast_generator(SPM_path, 'mvt', '')