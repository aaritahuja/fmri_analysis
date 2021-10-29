function task_second_level(MRI_directory, subjids, model_details, contrasts)

for i = 1:length(contrasts)
    confolders = {'simulation_over_baseline' 'perception_over_baseline' 'counting_over_baseline' 'simulation_over_counting' 'perception_over_counting'};
    directory = fullfile(MRI_directory, 'task_second_level', model_details, confolders{i}, date);
    %directory = fullfile(MRI_directory, 'task_second_level/uncertainty_parametric_orth/simulation_over_counting_both', date);
    if exist(directory, 'dir')
        warning('Directory already exists, making another')
        directories = dir(strcat(directory, '*'));
        directory = strcat(directory, '_', num2str(length(directories) + 1));
    end
    mkdir(directory)

    matlabbatch{1}.spm.stats.factorial_design.dir = {directory};
    filenames = {};
    for s = 1:length(subjids)
        filenames{s} = strcat(MRI_directory, '/', num2str(subjids(s)), '/models/', model_details, '/con_000', num2str(i), '.nii');
        %filenames{s} = strcat(MRI_directory, '/', num2str(subjids(s)), '/models/button_model_with_residuals/con_000', num2str(i), '.nii');
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
    
    %spm_jobman('interactive',matlabbatch)
    spm_jobman('run', matlabbatch)
    
end