%MRI_directory is the directory where all the MRI data is stored
%('/shared/lab/projects/analysis/aarit/MRI_data')
%subjid is the subject number
%good_runs is an integer list, in run numbering, of task runs to analyze
%localizer_runs is an integer list of localizer runs at the beginning of
%the session (not the end) to ignore

function SPM_directory = task_design_matrix(MRI_directory, subjid, task_runs, varargin)

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
        subject_directory = fullfile(MRI_directory, subjid, strcat('sub-glenn', sprintf('%03d', session)), 'ses-01');
                
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
    
    SPM_path = fullfile(subject_directory, 'models/task');
    SPM_directory = SPM_path;
    
    if ~exist(SPM_path, 'dir')
        mkdir(SPM_path)
    end
    
    n_runs = length(task_runs);
    bold_directory = fullfile(subject_directory, 'bold');
    behavior_directory = fullfile(subject_directory, 'behavior/task');
    folder_number = '%03d';
    matlabbatch{1}.spm.stats.fmri_spec.dir = {SPM_path};
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
        mvt = [];
        file_tracker = 0;
        run_repeat = [];
        runs = dir(fullfile(bold_directory, '*0*'));
        runs = struct2cell(runs);
        runs = runs(1, :);
        
        for i = 1:length(runs)
            run_directory = fullfile(bold_directory, runs(i));
            
            if i == 1
                finfo = dir(fullfile(run_directory{1}, 'rp*.txt'));
                tmpmvt = load(fullfile(finfo.folder, finfo.name));
                mvt = [tmpmvt];
            end
            
            if length(mvt) > 200
                file_list = dir(fullfile(run_directory{1}, 'swraf*'));
                run_repeat = [run_repeat; str2double(repmat(runs(i), length(file_list), 1))];
                file_tracker = file_tracker + length(file_list);
            end
            
        end
        
        mvt = [mvt run_repeat];
    end

    %General opening settings
    for i = 1:n_runs
        behavior_file = dir(fullfile(behavior_directory, strcat('*0', num2str(task_runs(i)), '.dgz')));
        behavior_group = dg_read(fullfile(behavior_directory, behavior_file.name));
        variant_column = cellstr(behavior_group.variant);
        variant = variant_column{1};
        
        %if strcmp(variant, 'native')
        %    error('Native variant detected')
        %end
        
        run_directory = fullfile(bold_directory, sprintf(folder_number, task_runs(i)));
        
        if isnumeric(subjid)
            file_list = dir(fullfile(run_directory, 'swraf*'));
            %file_list = dir(fullfile(run_directory, 'wraf*'));
        elseif ischar(subjid)
            file_list = dir(fullfile(run_directory, 's3wrlf*'));
        end
        
        files_with_path = strcat(run_directory, '/', {file_list.name});
        files_reshaped = reshape(files_with_path, [], 1);
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).scans = files_reshaped;
        
        %If no preprocessing package is specified, or if it is specified as
        %spm, do this
        if exist('preprocessing_package', 'var') == 0 || strcmp(preprocessing_package, 'spm')
            if length(mvt) > 200
                this_run_movement_idx = find(mvt(:, 7) == task_runs(i));
                this_run_movement_params = mvt(this_run_movement_idx, 1:6);
                for m = 1:6
                    matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress(m).name = sprintf('mvt%d', m);
                    matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress(m).val = this_run_movement_params(:, m);
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
        
        condition_list = {'pre_resp', 'first_two'};
        first_two_onsets = behavior_group.pre_resp_onsets(1:2);
        first_two_durations = behavior_group.pre_resp_durations(1:2);
        remaining_onsets = behavior_group.pre_resp_onsets(3:end);
        remaining_durations = behavior_group.pre_resp_durations(3:end);
        too_long = behavior_group.rts(3:end) >= ((prctile(behavior_group.rts(3:end), 75)) + (1.5 * iqr(behavior_group.rts(3:end)))) & behavior_group.rts(3:end) > 3000;
        skip_trials = behavior_group.rts(3:end) <= 0;
        if sum(too_long) > 0
            pre_resp_onsets = remaining_onsets(not(too_long));
            pre_resp_durations = remaining_durations(not(too_long));
            too_long_onsets = remaining_onsets(too_long);
            too_long_durations = remaining_durations(too_long);
            condition_list{end+1} = 'too_long';
        else
            pre_resp_onsets = remaining_onsets;
            pre_resp_durations = remaining_durations;
        end
        
        if sum(skip_trials) > 0
            pre_resp_onsets = pre_resp_onsets(not(skip_trials));
            pre_resp_durations = pre_resp_durations(not(skip_trials));
            skip_duration = ((behavior_group.stimoff(2) - behavior_group.stimon(2))/1000);
            skip_trial_onsets = remaining_onsets(skip_trials);
            skip_trial_durations = repmat(skip_duration, 1, length(skip_trial_onsets));
            condition_list{end+1} = 'skip_trials';
        end
        
        %Condition 1: pre-response period
        if strcmp(variant, 'native')
            condition_name = strcat(variant, '_pre_response');
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(1).name = condition_name;
        else
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(1).name = variant;
        end

        condition_number = find(strcmp(condition_list, 'pre_resp'));
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).onset = pre_resp_onsets;
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).duration = pre_resp_durations;
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).orth = 1;
        
        %Condition 2: First two trials to skip
        condition_number = find(strcmp(condition_list, 'first_two'));
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).name = 'first_trials';
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).onset = first_two_onsets;
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).duration = first_two_durations;
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).orth = 1;
        
        if sum(too_long) > 0
            %Condition n: Trials with outlier reaction times to skip
            condition_number = find(strcmp(condition_list, 'too_long'));
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).name = 'long_trials';
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).onset = too_long_onsets;
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).duration = too_long_durations;
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).orth = 1;
        end
        
        if sum(skip_trials) > 0
            %Condition n: Trials that were skipped
            condition_number = find(strcmp(condition_list, 'skip_trials'));
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).name = 'skip_trials';
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).onset = skip_trial_onsets;
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).duration = skip_trial_durations;
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).orth = 1;
        end
        
        %Condition n: Post-response period, only if in native variant
        if strcmp(variant, 'native')
            condition_list{end+1} = 'post_resp';
            condition_number = find(strcmp(condition_list, 'post_resp'));
            condition_name = strcat(variant, '_post_response');
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).name = condition_name;
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).onset = behavior_group.response_onsets;
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).duration = behavior_group.post_resp_duration;
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(condition_number).orth = 1;
        end

        
        %General closing settings
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi_reg = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).hpf = 128;
    end
%%        
%Uncertainty parametric regressor
%         if include_uncertainty == 1
%             if strcmp(variant, 'simulation')
%                 matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(1).pmod.name = 'uncertainty';
%                 matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(1).pmod.poly = 1;
%                 matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(1).pmod.param = behavior_group.uncertainty;
%             else
%                matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
%             end
%         end
%%      
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    
    if isnumeric(subjid)
        matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    elseif ischar(subjid)
        if exist('preprocessing_package', 'var') == 0 || strcmp(preprocessing_package, 'spm')
            if exist(fullfile(subject_directory, 'anatomy/wssbrainmask.nii'), 'file')
                matlabbatch{1,1}.spm.stats.fmri_spec.mask = ...
                    {fullfile(subject_directory, 'anatomy/wssbrainmask.nii')};
            else
                error('ms_1stLevel:mask', ...
                    'Subject %d does not have wssbrainmask.nii in anatomy folder', sessnum)
            end
        elseif exist('preprocessing_package', 'var') && strcmp(preprocessing_package, 'afni')
            if exist(fullfile(subject_directory, 'warped/NMT_v2.0_sym_05mm_brainmask.nii'), 'file')
                matlabbatch{1,1}.spm.stats.fmri_spec.mask = ...
                    {fullfile(subject_directory, 'warped/NMT_v2.0_sym_05mm_brainmask.nii')};
            else
                error('ms_1stLevel:mask', ...
                    'Subject %d does not have NMT_v2.0_sym_05mm_brainmask.nii in warped folder - might need to be expanded from .gz', sessnum)
            end
        end
    end

    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    %matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    %matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    %matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    spm_jobman('interactive',matlabbatch)
    %spm_jobman('run', matlabbatch)
end