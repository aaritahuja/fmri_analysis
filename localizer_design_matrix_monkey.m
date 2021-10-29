function SPM_path = localizer_design_matrix_monkey(MRI_directory, subjid, run_table, hrf_type)
%%
    sessions = transpose(unique(run_table.session));
    for i = 1:length(sessions)
        session_runs = transpose(run_table.run(run_table.session == sessions(i)));
        localizer_runs{i} = session_runs;
    end
    
    %Initialize SPM
    spm fmri
    spm_jobman('initcfg');
    clear matlabbatch
    %addpath('/shared/lab/projects/analysis/aarit/MATLAB')
    conditions = ["motion", "flicker", "static"];
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
    SPM_path = fullfile(MRI_directory, subjid, 'localizer_models/spm', combined_session);
    if ~exist(SPM_path, 'dir')
        mkdir(SPM_path)
    end
    
    %copyfile(fullfile(MRI_directory, 'Glenn/task_models/NMT_v2.0_sym_05mm_brainmask.nii'), SPM_path)
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
        error('TR unclear. Check scan sequence')
    end
    
    overall_run = 0; 
    for s = 1:length(sessions)
        session = sessions(s);
        subject_directory = fullfile(MRI_directory, subjid, strcat('sub-glenn', sprintf('%03d', session)), 'ses-01');
        localizer_run = localizer_runs{s};
        n_runs = length(localizer_run);
        bold_directory = fullfile(subject_directory, 'bold');
        behavior_directory = fullfile(subject_directory, 'behavior/localizer');
        folder_number = '%03d';
        
        for i = 1:n_runs
            overall_run = overall_run + 1;
            behavior_file = dir(fullfile(behavior_directory, strcat('*0', num2str(localizer_run(i)), '.dgz')));
            behavior_group = dg_read(fullfile(behavior_directory, behavior_file.name));
            
            %Catch for cases where the list of obs_times might be too long, for
            %e.g if a button was pressed during the localizer block
            if length(behavior_group.obs_times) ~= length(behavior_group.id)
                disp('mismatch detected')
                obs_diff = diff(behavior_group.obs_times);
                short_indices = find(obs_diff < 18000);
                long_indices = find(obs_diff > 24000);
                bad_indices = [short_indices  long_indices];
                behavior_group.obs_times(bad_indices) = [];
                %behavior_group.obs_times = behavior_group.obs_times(short_indices);
                if length(behavior_group.obs_times) ~= length(behavior_group.id)
                    error('Unable to resolve discrepancy between list of obs_times and list of ids, fix manually');
                end
            end
            
            onset_times = (behavior_group.obs_times)/1000;
            variant_ids = behavior_group.variant_id;
            
            run_directory = fullfile(bold_directory, sprintf(folder_number, localizer_run(i)));
            file_list = dir(fullfile(run_directory, 's3wrlf*'));
            
            files_with_path = strcat(run_directory, '/', {file_list.name});
            files_reshaped = reshape(files_with_path, [], 1);
            matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).scans = files_reshaped;
            
            finfo = dir(fullfile(bold_directory, sprintf(folder_number, localizer_run(i)), 'mot_demean*'));
            mvt = load(fullfile(finfo.folder, finfo.name));
            for m = 1:6
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).regress(m).name = sprintf('mvt%d', m);
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).regress(m).val = nonzeros(mvt(:, m));
            end
            
            %I made a change to how the variant ids are coded in the first week of
            %December 2019, so this lets them both coexist
            run_date = datetime(behavior_file.date);
            if run_date < datetime('19-Dec-2019')
                variant_id_list = [1 2 3];
            else
                variant_id_list = [5 6 7];
            end
            
            %Set the three conditions
            for c = 1:length(conditions)
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(c).name = conditions{c};
                logical_selector = (variant_ids == variant_id_list(c));
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(c).onset = onset_times(logical_selector);
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(c).duration = 20;
                
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(c).tmod = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).cond(c).orth = 1;
            end

            %General closing settings
            matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).multi = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).multi_reg = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess(overall_run).hpf = 128;
        end
    end

    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.0;
    
 
    matlabbatch{1,1}.spm.stats.fmri_spec.mask = {fullfile(SPM_path, 'resampled_mask.nii')};
   
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    %spm_jobman('interactive',matlabbatch)
    spm_jobman('run', matlabbatch)
end