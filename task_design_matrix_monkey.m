%MRI_directory is the directory where all the MRI data is stored
%('/shared/lab/projects/analysis/aarit/MRI_data')
%subjid is the subject number
%good_runs is an integer list, in run numbering, of task runs to analyze
%localizer_runs is an integer list of localizer runs at the beginning of
%the session (not the end) to ignore

function SPM_path = task_design_matrix_monkey(MRI_directory, subjid, run_table, hrf_type, censor_bad_trials)

    sessions = transpose(unique(run_table.session));
    for i = 1:length(sessions)
        session_runs = transpose(run_table.run(run_table.session == sessions(i)));
        task_runs{i} = session_runs;
    end
    
    %Initialize SPM
    spm fmri
    spm_jobman('initcfg');
    clear matlabbatch
    
    if ismac
        hrf_path = fullfile('/Users/Aaru/Documents/MATLAB', hrf_type);
    elseif isunix
        hrf_path = fullfile('/home/aarit/MATLAB', hrf_type);
    else
        hrf_path = fullfile('C:/Users/lab/Documents/MATLAB', hrf_type);
    end
    
    %Adding the hrf_path to the top will ensure that spm picks the appropriate scripts
    %from this folder
    addpath(hrf_path, '-begin')
    
    combined_session = string(sessions);
    combined_session = [combined_session{:}];
    SPM_path = fullfile(MRI_directory, subjid, 'task_models/spm', combined_session);
    if ~exist(SPM_path, 'dir')
        mkdir(SPM_path)
    end
    copyfile(fullfile(MRI_directory, 'Glenn/sub-glenn005/ses-01/warped/NMT_v2.0_sym_05mm_brainmask.nii'), SPM_path)
    copyfile(fullfile(MRI_directory, 'Glenn/task_models/resampled_mask.nii'), SPM_path)
    
    matlabbatch{1}.spm.stats.fmri_spec.dir = {SPM_path};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    
    if hrf_type == 'hrf_mion'
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.8;
    elseif hrf_type == 'hrf_bold'
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2.76;
    else
        disp('TR unclear. Check scan sequence')
    end
   
   overall_run = 0; 
   for s = 1:length(sessions)
       session = sessions(s);
       subject_directory = fullfile(MRI_directory, subjid, strcat('sub-glenn', sprintf('%03d', session)), 'ses-01');
       task_run = task_runs{s};
       n_runs = length(task_run);
       bold_directory = fullfile(subject_directory, 'bold_scaled');
       behavior_directory = fullfile(subject_directory, 'behavior/task');
       folder_number = '%03d';
       
    %General opening settings
    for i = 1:n_runs
        overall_run = overall_run + 1;
        behavior_file = dir(fullfile(behavior_directory, strcat('*0', num2str(task_run(i)), '.dgz')));
        behavior_group = dg_read(fullfile(behavior_directory, behavior_file.name));
        variant_column = cellstr(behavior_group.variant);
        variant = variant_column{1};
        
        
        run_directory = fullfile(bold_directory, sprintf(folder_number, task_run(i)));
        file_list = dir(fullfile(run_directory, 's3wrlf*'));
        
        files_with_path = strcat(run_directory, '/', {file_list.name});
        files_reshaped = reshape(files_with_path, [], 1);
        matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).scans = files_reshaped;
        
        finfo = dir(fullfile(bold_directory, sprintf(folder_number, task_run(i)), 'mot_demean*'));
        mvt = load(fullfile(finfo.folder, finfo.name));
        
        for m = 1:6
            matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).regress(m).name = sprintf('mvt%d', m);
            matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).regress(m).val = nonzeros(mvt(:, m));
        end
%%        
%         censor_finfo = dir(fullfile(subject_directory, 'preprocessed/motion_*censor.1D'));
%         censor = load(fullfile(censor_finfo.folder, censor_finfo.name));
%         run_censor = not(censor(run_indices))*-1;
%         cluster = [];
%         n_clusters = 0;
%         for c = 1:length(run_censor)
%             if run_censor(c) == 0
%                 cluster(c) = 0;
%             elseif run_censor(c) == -1
%                 if run_censor(c-1) ~= -1
%                     n_clusters = n_clusters + 1;
%                 end
%                 cluster(c) = n_clusters;
%             end
%         end
%         cluster = transpose(cluster);
%         
%         
%         for c = 1:n_clusters
%             scrub_regressor = (cluster == c) * -1;
%             matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).regress(m+c).name = sprintf('censor%d', c);
%             matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).regress(m+c).val = scrub_regressor;
%         end
% 
%%
        run_indices = find(mvt(:, 1) ~= 0);
        
        preprocessed_directory = fullfile(subject_directory, 'preprocessed');
        censor_finfo = dir(fullfile(preprocessed_directory, 'motion_*censor.1D'));
        censor = load(fullfile(censor_finfo.folder, censor_finfo.name));
        run_censor = not(censor(run_indices));
        run_censor = double(run_censor);
%%        
        if censor_bad_trials  == 1
            %Get scan indices for scans spanning the start of the first trial to the start of the third trial
            first_two_trials = behavior_group.pre_resp_onsets([1, 3]);
            first_scan_times = first_two_trials/1.8;
            first_scan_indices = transpose([round(first_scan_times(1)):floor(first_scan_times(2) - 1)]);
        
            too_long_scan_indices = [];
            skip_trial_scan_indices = [];
            %Get scan indices for trials with outlier reaction times
            too_long = find(behavior_group.rts(3:end) >= ((prctile(behavior_group.rts(3:end), 75)) + (1.5 * iqr(behavior_group.rts(3:end)))) & behavior_group.rts(3:end) > 4000);
            too_long = too_long + 2;
            for tl = 1:length(too_long)
                too_long_trial_start = behavior_group.pre_resp_onsets(too_long(tl));
                too_long_trial_end = behavior_group.pre_resp_onsets(too_long(tl)) + (behavior_group.stimoff(too_long(tl))/1000) + (behavior_group.feedback_iti(too_long(tl))/1000);
                too_long_scans = transpose([round((too_long_trial_start)/1.8):floor((too_long_trial_end)/1.8)]);
                too_long_scan_indices = [too_long_scan_indices; too_long_scans];
            end
        
            %Get scan indices for skip trials
            skip_trials = find(behavior_group.rts > 12000);
            for st = 1:length(skip_trials)
                skip_trial_start = behavior_group.pre_resp_onsets(skip_trials(st));
                skip_trial_end = behavior_group.pre_resp_onsets(skip_trials(st)) + (behavior_group.stimoff(skip_trials(st))/1000);
                skip_trial_scans = transpose([round((skip_trial_start)/1.8):round((skip_trial_end)/1.8)]);
                skip_trial_scan_indices = [skip_trial_scan_indices; skip_trial_scans];
            end
        
            %Manually add in censoring for the first to trials, skip trials, and
            %outlier rt trials. Then concat this array to cumulative array
            run_censor([skip_trial_scan_indices; too_long_scan_indices; first_scan_indices]) = 1;
            
            condition_list = {'pre_resp'};
            
            pre_resp_onsets = behavior_group.pre_resp_onsets;
            pre_resp_durations = behavior_group.pre_resp_durations;
%%          
        elseif censor_bad_trials == 0
            
            condition_list = {'pre_resp', 'first_two'};
            first_two_onsets = behavior_group.pre_resp_onsets(1:2);
            first_two_durations = behavior_group.pre_resp_durations(1:2);
            
            remaining_onsets = behavior_group.pre_resp_onsets(3:end);
            remaining_durations = behavior_group.pre_resp_durations(3:end);
            remaining_resp = behavior_group.resp(3:end);
            remaining_rts = behavior_group.rts(3:end);
            
            skip_trials = remaining_resp == -1;
            rts_to_check = remaining_rts(remaining_resp > -1);
            too_long = rts_to_check >= ((prctile(rts_to_check, 75)) + (1.5 * iqr(rts_to_check))) & rts_to_check > 4000;
            
            if sum(too_long) > 0
                long_rts = rts_to_check(find(too_long));
                [tf,idx] = ismember(remaining_rts, long_rts);
                
                pre_resp_onsets = remaining_onsets(not(tf));
                pre_resp_durations = remaining_durations(not(tf));
                too_long_onsets = remaining_onsets(tf);
                too_long_durations = remaining_durations(tf);
                condition_list{end+1} = 'too_long';
            else
                pre_resp_onsets = remaining_onsets;
                pre_resp_durations = remaining_durations;
            end
            
            if sum(skip_trials) > 0
                skip_onsets = remaining_onsets(find(skip_trials));
                [tf,idx] = ismember(pre_resp_onsets, skip_onsets);
                
                skip_trial_onsets = pre_resp_onsets(tf);
                skip_trial_durations = pre_resp_durations(tf);
                
                pre_resp_onsets = pre_resp_onsets(not(tf));
                pre_resp_durations = pre_resp_durations(not(tf));
                condition_list{end+1} = 'skip_trials';
            end
            
            combine_error_trials = 0;
            if combine_error_trials == 1
                error_onsets = first_two_onsets;
                error_durations = first_two_durations;
                
                if sum(too_long) > 0
                    error_onsets = [error_onsets; too_long_onsets];
                    error_durations = [error_durations; too_long_durations];
                end
                
                if sum(skip_trials) > 0
                    error_onsets = [error_onsets; skip_trial_onsets];
                    error_durations = [error_durations; skip_trial_durations];
                end
                
                condition_list = {'pre_resp', 'error'};
            end
            
            %Condition 1: pre-response period
            if strcmp(variant, 'native')
                condition_name = strcat(variant, '_pre_response');
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(1).name = condition_name;
            else
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(1).name = variant;
            end
            
        else
            error('Unclear what to do with bad trials')
        end
        
        censor_indices = find(run_censor == 1);
        for c = 1:length(censor_indices)
            this_censor = zeros(length(run_censor), 1);
            this_censor(censor_indices(c)) = 1;
            matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).regress(m+c).name = sprintf('censor%d', c);
            matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).regress(m+c).val = this_censor;
        end
%%      

        condition_number = find(strcmp(condition_list, 'pre_resp'));
        matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).name = variant;
        matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).onset = pre_resp_onsets;
        matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).duration = pre_resp_durations;
        matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).orth = 1;
%%      
        if censor_bad_trials == 0
            
            if combine_error_trials == 0
                %Condition 2: First two trials to skip
                condition_number = find(strcmp(condition_list, 'first_two'));
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).name = 'first_trials';
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).onset = first_two_onsets;
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).duration = first_two_durations;
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).tmod = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).pmod = struct('name', {}, 'param', {}, 'poly', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).orth = 1;
                
                if sum(too_long) > 0
                    %Condition n: Trials with outlier reaction times to skip
                    condition_number = find(strcmp(condition_list, 'too_long'));
                    matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).name = 'long_trials';
                    matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).onset = too_long_onsets;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).duration = too_long_durations;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).tmod = 0;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).pmod = struct('name', {}, 'param', {}, 'poly', {});
                    matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).orth = 1;
                end
                
                if sum(skip_trials) > 0
                    %Condition n: Trials that were skipped
                    condition_number = find(strcmp(condition_list, 'skip_trials'));
                    matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).name = 'skip_trials';
                    matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).onset = skip_trial_onsets;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).duration = skip_trial_durations;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).tmod = 0;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).pmod = struct('name', {}, 'param', {}, 'poly', {});
                    matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).orth = 1;
                end
            elseif combine_error_trials == 1
                %Condition 2: Error trials
                condition_number = find(strcmp(condition_list, 'error'));
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).name = 'error_trials';
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).onset = error_onsets;
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).duration = error_durations;
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).tmod = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).pmod = struct('name', {}, 'param', {}, 'poly', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).orth = 1;
            end
            
            %Condition n: Post-response period, only if in native variant
            if strcmp(variant, 'native')
                condition_list{end+1} = 'post_resp';
                condition_number = find(strcmp(condition_list, 'post_resp'));
                condition_name = strcat(variant, '_post_response');
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).name = condition_name;
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).onset = behavior_group.response_onsets;
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).duration = behavior_group.post_resp_duration;
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).tmod = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).pmod = struct('name', {}, 'param', {}, 'poly', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(condition_number).orth = 1;
            end
        end
%%        
        %General closing settings
        matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).multi = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).multi_reg = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).hpf = 128;
    end
   end
%%        
  
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0;

    matlabbatch{1,1}.spm.stats.fmri_spec.mask = {fullfile(SPM_path, 'resampled_mask.nii')};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    %spm_jobman('interactive',matlabbatch)
    spm_jobman('run', matlabbatch)
end