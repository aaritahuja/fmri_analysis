function full_table = run_table_localizer(MRI_directory, sessions)

session = [];
run = [];
percent_censor = [];
usable = [];


for i = 1:length(sessions)
    sess = sprintf('%03d', sessions(i));
    session_directory = fullfile(MRI_directory, strcat('Glenn/sub-glenn', sess), 'ses-01');
    bold_directory = strcat(session_directory, '/bold/');
    behavior_directory = strcat(session_directory, '/behavior/localizer/');
    bold_runs = struct2cell(dir(strcat(bold_directory, '0*')));
    runs = str2double(bold_runs(1, :));
    
    for r = 1:length(runs)
        behavior_file = dir(fullfile(behavior_directory, strcat('*motionloc*', sprintf('%02d', runs(r)), '.dgz')));
        if length(behavior_file) > 0
            finfo = dir(fullfile(bold_directory, sprintf('%03d', runs(r)), 'mot_demean*'));
            mvt = load(fullfile(finfo.folder, finfo.name));
            run_indices = find(mvt(:, 1) ~= 0);
            
            censor_finfo = dir(fullfile(session_directory, 'preprocessed/motion_*censor.1D'));
            censor = load(fullfile(censor_finfo.folder, censor_finfo.name));
            run_censor = not(censor(run_indices));
            proportion = (sum(run_censor)/length(run_censor))*100;
            
            behavior_group = dg_read(fullfile(behavior_file.folder, behavior_file.name));
                        
            if proportion <= 15
                use = 'Y';
            else
                use = 'N';
            end
            
            session = [session; sessions(i)];
            run = [run; runs(r)];
            percent_censor = [percent_censor; proportion];
            usable = [usable; use];
        end
    end
end
      
full_table = table(session, run, percent_censor, usable);

        
        
    