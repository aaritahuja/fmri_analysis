function localizer_second_level(MRI_directory, subjids, contrasts, varargin)

for i = 1:length(contrasts)
    contrast = contrasts(i);
    if contrast == 1
        conname = 'motion_over_static';
    elseif contrast == 2
        conname = 'static_over_motion';
    elseif contrast == 3
        conname = 'flicker_over_static';
    end
  
    if isnumeric(subjids(1))
        directory = fullfile(MRI_directory, 'localizer_second_level', conname, date);
        if exist(directory, 'dir')
            warning('Directory already exists, making another')
            directories = dir(strcat(directory, '*'));
            directory = strcat(directory, '_', num2str(length(directories) + 1));
        end
        mkdir(directory)
    elseif ischar(subjids)
        sessions = varargin{1};
        directory = fullfile(MRI_directory, subjids, 'localizer_second_level', conname, date);
        if exist(directory, 'dir')
            warning('Directory already exists, making another')
            directories = dir(strcat(directory, '*'));
            directory = strcat(directory, '_', num2str(length(directories) + 1));
        end
        mkdir(directory)
    end

    matlabbatch{1}.spm.stats.factorial_design.dir = {directory};
    filenames = {};
    if isnumeric(subjids(1))
        for s = 1:length(subjids)
            filenames{s} = strcat(MRI_directory, '/', num2str(subjids(s)), '/models/localizer/con_000', num2str(contrast), '.nii');
            
        end
    elseif ischar(subjids)
        for s = 1:length(sessions)
            filenames{s} = fullfile(MRI_directory, subjids, sessions{s}, strcat('ses-01/models/localizer/con_000', num2str(contrast), '.nii'));
        end
    end
    
    filenames = reshape(filenames, [length(filenames), 1]);
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = filenames;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    %Estimation
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 1;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    spm_jobman('interactive',matlabbatch)
    %spm_jobman('run', matlabbatch)
end