function second_level_contrasts (MRI_directory, model_details, contrasts)

for i = contrasts
    confolders = {'simulation_over_baseline' 'perception_over_baseline' 'counting_over_baseline' 'simulation_over_counting' 'perception_over_counting'};  
    directory = fullfile(MRI_directory, 'task_second_level', model_details, confolders{i}, date);
    %make the contrast
    matlabbatch{1}.spm.stats.con.spmmat = {fullfile(directory, 'SPM.mat')};
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = confolders{i} ;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1];
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.delete = 0;
    %spm_jobman('interactive',matlabbatch)
    spm_jobman('run',matlabbatch)
end
