function condition_movement_regressors(MRI_directory, subjid, task_runs, contrasts)

subject_directory = fullfile(MRI_directory, num2str(subjid));
behavior_directory = fullfile(subject_directory, 'behavior/task');
SPM_path = fullfile('/Users/Aaru/Documents/MRI_data', num2str(subjid), 'models/task_with_duration');
load(fullfile(SPM_path, 'SPM.mat'));

for i = 1:length(task_runs)
    behavior_file = dir(fullfile(behavior_directory, strcat('*0', num2str(task_runs(i)), '.dgz')));
    behavior_group = dg_read(fullfile(behavior_directory, behavior_file.name));
    variant_column = cellstr(behavior_group.variant);
    variant = variant_column{1};
    variants{i} = variant;
end

for c = 1:length(contrasts)
    contrast_blocks = find(double(contains(variants, contrasts{c})));    
    con_mvt = zeros(1, length(SPM.xX.name));
    
   for i = contrast_blocks
        con_m = double(contains(SPM.xX.name, strcat('Sn(', num2str(i), ') mvt')));
        con_mvt = con_mvt + con_m;
   end
   
   contrast_name = strcat(contrasts{c}, '_mvt');
   matlabbatch{1}.spm.stats.con.spmmat = {fullfile(SPM_path, 'SPM.mat')};
   matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = contrast_name ;
   matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = con_mvt;
   matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
   matlabbatch{1}.spm.stats.con.delete = 0;
   %spm_jobman('interactive',matlabbatch)
   spm_jobman('run',matlabbatch)
    
end
  
end