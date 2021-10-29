function SPM_directory = localizer_design_matrix(MRI_directory, subjid, localizer_runs, varargin)
%%
    %Initialize SPM
    spm fmri
    spm_jobman('initcfg');
    clear matlabbatch
    %addpath('/shared/lab/projects/analysis/aarit/MATLAB')
%%
    n_runs = length(localizer_runs);
    conditions = ["motion", "flicker", "static"];
    
    if isnumeric(subjid) %If subjid is a number, we assume human
        subject_directory = fullfile(MRI_directory, num2str(subjid));
    elseif ischar(subjid) %Else, if subjid is a string, we assume monkey
        %Grab session number from varagin
        if nargin > 0
            session = varargin{1};
        else
            error('Session number not specified')
        end
        
        %Grab hrf type from varargin, if specified
        if nargin > 1
            hrf_type = varargin{2};
            if ismac
                hrf_path = fullfile('/Users/Aaru/Documents/MATLAB', hrf_type);
            elseif isunix
                hrf_path = fullfile('/home/aarit/MATLAB', hrf_type);
            else
                hrf_path = fullfile('C:/Users/lab/Documents/MATLAB', hrf_type);
            end
        end
        
        subject_directory = fullfile(MRI_directory, subjid, session);
        
        %Adding the hrf_path to the top will ensure that spm picks the appropriate scripts
        %from this folder
        addpath(hrf_path, '-begin')
        
        if length(varargin) == 3
           preprocessing_package = varargin{3};
        end
        
    end
    
    SPM_path = fullfile(subject_directory, 'models/localizer');
    SPM_directory = SPM_path;
    
    if ~exist(SPM_path, 'dir')
        mkdir(SPM_path)
    end
    
    bold_directory = fullfile(subject_directory, 'bold');
    behavior_directory = fullfile(subject_directory, 'behavior/localizer');
    folder_number = '%03d';
    matlabbatch{1}.spm.stats.fmri_spec.dir = {SPM_path};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    
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
    
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
 
    if exist('preprocessing_package', 'var') == 0 || strcmp(preprocessing_package, 'spm')
        %Value to update to track cumulative files after each run
        file_tracker = 0;
        
        %Get movement regressors (the current format assumes realignment was done all at once, which means there is only one file for all runs
        mvt = [];
        finfo = dir(fullfile(bold_directory, sprintf(folder_number, localizer_runs(1)), 'rp*.txt'));
        if ~isempty(finfo)
            %tmpmvt = load(fullfile(bold_directory, sprintf(folder_number, runs(1)), finfo.name));
            tmpmvt = load(fullfile(finfo.folder, finfo.name));
            mvt = [mvt; tmpmvt];
        else
            disp('Problem getting movement regressors. Maybe the realignment for all runs was not done at the same time?')
        end
    end
%%
    for i = 1:n_runs
        %General opening settings
        behavior_file = dir(fullfile(behavior_directory, strcat('*0', num2str(localizer_runs(i)), '.dgz')));
        behavior_group = dg_read(fullfile(behavior_directory, behavior_file.name));
        
        %Catch for cases where the list of obs_times might be too long, for
        %e.g if a button was pressed during the localizer block
        if length(behavior_group.obs_times) ~= length(behavior_group.id)
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
        run_directory = fullfile(bold_directory, sprintf(folder_number, localizer_runs(i)));
        
        %Different file prefix for monkey vs human
        if isnumeric(subjid)
            file_list = dir(fullfile(run_directory, 'swraf*'));
            if isempty(file_list)
                file_list = dir(fullfile(run_directory, 'wraf*'));
            end
        elseif ischar(subjid)
            file_list = dir(fullfile(run_directory, 's3wrlf*'));
        end
        
        files_with_path = strcat(run_directory, '/', {file_list.name});
        files_reshaped = reshape(files_with_path, [], 1);
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).scans = files_reshaped;
 
        %If no preprocessing package is specified, or if it is specified as
        %spm, do this
        if exist('preprocessing_package', 'var') == 0 || strcmp(preprocessing_package, 'spm')
            % Quick and dirty check for whether realignment was done all at once or run by run
            if length(mvt) > 200
                %error('Have not yet figured out how to grab movement regressors for localizer when task runs are in between')
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
                finfo = dir(fullfile(bold_directory, sprintf(folder_number, good_runs(i)), 'rp*.txt'));
                mvt = load(fullfile(finfo.folder, finfo.name));
                %Add movmement regressors 6 times, once for each column in mvt which correspond to one of 6 motion parameters
                for m = 1:6
                    matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress(m).name = sprintf('mvt%d', m);
                    matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress(m).val = mvt(:, m);
                end
            end
        elseif exist('preprocessing_package', 'var') && strcmp(preprocessing_package, 'afni')
            finfo = dir(fullfile(bold_directory, sprintf(folder_number, localizer_runs(i)), 'mot_demean*'));
            mvt = load(fullfile(finfo.folder, finfo.name));
            for m = 1:6
                matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress(m).name = sprintf('mvt%d', m);
                matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress(m).val = nonzeros(mvt(:, m));
            end
        end

%%
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
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(c).name = conditions{c};
            logical_selector = (variant_ids == variant_id_list(c));
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(c).onset = onset_times(logical_selector);            
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(c).duration = 20;

            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(c).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(c).orth = 1;
        end
%%      
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
    
    %Different mask settings for monkey vs human
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
                    'Subject does not have NMT_v2.0_sym_05mm_brainmask.nii in warped folder - might need to be expanded from .gz')
            end
        end
    end
    
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    %spm_jobman('interactive',matlabbatch)
    spm_jobman('run', matlabbatch)
end