%MRI_directory is the directory where all the MRI data is stored
%('/shared/lab/projects/analysis/aarit/MRI_data')
%subjid is the subject number
%good_runs is an integer list, in run numbering, of task runs to analyze
%localizer_runs is an integer list of localizer runs at the beginning of
%the session (not the end) to ignore

function SPM_directory = button_design_matrix(MRI_directory, subjid, task_runs, localizer_runs, modeltype, varargin)

    %Initialize SPM
    spm fmri
    spm_jobman('initcfg');
    clear matlabbatch

    if isnumeric(subjid)
        subject_directory = fullfile(MRI_directory, num2str(subjid));
    elseif ischar(subjid)
        if length(varargin) > 0
            session = varargin{1};
        else
            error('Session number not specified')
        end
        subject_directory = fullfile(MRI_directory, subjid, session);
                
        %Grab hrf type from varargin, if specified
        if length(varargin) < 2
            error('Subject is monkey, but HRF type is not specified')
        else
            hrf_type = varargin{2};
            if ismac
                hrf_path = fullfile('/Users/Aaru/Documents/MATLAB', hrf_type);
            elseif isunix
                hrf_path = fullfile('/home/aarit/MATLAB', hrf_type);
            else
                hrf_path = fullfile('C:/Users/lab/Documents/MATLAB', hrf_type);
            end
        end
        %Adding the hrf_path to the top will ensure that spm picks the appropriate scripts
        %from this folder
        addpath(hrf_path, '-begin')
       if length(varargin) == 3
           preprocessing_package = varargin{3};
       end
    end
    
   
    n_runs = length(task_runs);
    bold_directory = fullfile(subject_directory, 'bold');
    behavior_directory = fullfile(subject_directory, 'behavior/task');
    folder_number = '%03d';
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    
    %Different TRs for human and monkey
    if isnumeric(subjid)
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
    elseif ischar(subjid)
        if hrf_type == 'hrf_mion'
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.8;
        elseif hrf_type == 'hrf_bold'
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2.76;
        else
            disp('TR unclear. Check scan sequence')
        end
    end
   
    %If no preprocessing package is specified, or if it is specified as
    %spm, do this because spm makes one large file
    if exist('preprocessing_package', 'var') == 0 || strcmp(preprocessing_package, 'spm')
        %Value to update to track cumulative files after each run
        file_tracker = 0;
        
        for i = 1:length(localizer_runs)
            localizer_directory = fullfile(bold_directory, sprintf(folder_number, localizer_runs(i)));
            if isnumeric(subjid)
                localizer_files = dir(fullfile(localizer_directory, 'swraf*'));
            elseif ischar(subjid)
                localizer_files = dir(fullfile(localizer_directory, 's3wrlf*'));
            end
            file_tracker = file_tracker + length(localizer_files);
        end
        
        %Get spm movement regressors (the current format assumes realignment was done all at once, which means there is only one file for all runs
        mvt = [];
        finfo = dir(fullfile(bold_directory, sprintf(folder_number, localizer_runs(1)), 'rp*.txt'));
        if ~isempty(finfo)
            tmpmvt = load(fullfile(finfo.folder, finfo.name));
            mvt = [mvt; tmpmvt];
        else
            disp('Problem getting movement regressors. No realignment text file in first localizer run')
        end
    end

    %General opening settings
    for i = 1:n_runs
        behavior_file = dir(fullfile(behavior_directory, strcat('*0', num2str(task_runs(i)), '.dgz')));
        behavior_group = dg_read(fullfile(behavior_directory, behavior_file.name));
        variant_column = cellstr(behavior_group.variant);
        variant = variant_column{1};
        run_directory = fullfile(bold_directory, sprintf(folder_number, task_runs(i)));
        if isnumeric(subjid)
            file_list = dir(fullfile(run_directory, 'wraf*'));
        elseif ischar(subjid)
            file_list = dir(fullfile(run_directory, 'wrlf*'));
        end
        files_with_path = strcat(run_directory, '/', {file_list.name});
        files_reshaped = reshape(files_with_path, [], 1);
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).scans = files_reshaped;
        
        %If no preprocessing package is specified, or if it is specified as
        %spm, do this
        if exist('preprocessing_package', 'var') == 0 || strcmp(preprocessing_package, 'spm')
            % Quick and dirty check for whether realignment was done all at once or run by run
            if length(mvt) > 200
                %Add movmement regressors 6 times, once for each column in mvt which correspond to one of 6 motion parameters
                mov_idx = (file_tracker + 1): (file_tracker + length(files_with_path));
                file_tracker = file_tracker + length(files_with_path);
                %disp(file_tracker);
                for m = 1:6
                    matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress(m).name = sprintf('mvt%d', m);
                    matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress(m).val = mvt(mov_idx, m);
                end
            else
                disp('Length of movement regressors in first file is very short. Attempting to get run by run movement regressors')
                finfo = dir(fullfile(bold_directory, sprintf(folder_number, task_runs(i)), 'rp*.txt'));
                mvt = load(fullfile(finfo.folder, finfo.name));
                %Add movmement regressors 6 times, once for each column in mvt which correspond to one of 6 motion parameters
                for m = 1:6
                    matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress(m).name = sprintf('mvt%d', m);
                    matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress(m).val = mvt(:, m);
                end
            end
        elseif exist('preprocessing_package', 'var') && strcmp(preprocessing_package, 'afni')
            finfo = dir(fullfile(bold_directory, sprintf(folder_number, task_runs(i)), 'mot_demean*'));
            mvt = load(fullfile(finfo.folder, finfo.name));
            for m = 1:6
                matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress(m).name = sprintf('mvt%d', m);
                matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress(m).val = nonzeros(mvt(:, m));
            end
        end

        if strcmp(modeltype, 'ROI')
            SPM_path = fullfile(subject_directory, 'models/button_ROI_model');
            condition_name = 'button_press';
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(1).name = condition_name;
        elseif strcmp(modeltype, 'variant_specific')
            SPM_path = fullfile(subject_directory, 'models/button_variant_model');
            condition_name = strcat(variant, '_button');
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(1).name = condition_name;
        else
            error('Type of model not specified - is this variant specific or just a general ROI?')
        end
        
        SPM_directory = SPM_path;
        if ~exist(SPM_path, 'dir')
            mkdir(SPM_path)
        end
        matlabbatch{1}.spm.stats.fmri_spec.dir = {SPM_path};
        
        onsets = behavior_group.response_onsets;
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(1).onset = onsets;
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(1).duration = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(1).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(1).orth = 1;
        
        %General closing settings
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi_reg = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).hpf = 128;
    end

    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 1;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    %spm_jobman('interactive',matlabbatch)
    spm_jobman('run', matlabbatch)
end