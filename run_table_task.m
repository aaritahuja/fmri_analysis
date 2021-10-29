function [full_table, usable_table] = run_table(MRI_directory, sessions)

session = [];
run = [];
variant = [];
percent_correct = [];
side_bias = [];
n_skip_trials = [];
n_too_long = [];
percent_censor = [];
usable = [];


for i = 1:length(sessions)
    sess = sprintf('%03d', sessions(i));
    session_directory = fullfile(MRI_directory, strcat('Glenn/sub-glenn', sess), 'ses-01');
    bold_directory = strcat(session_directory, '/bold/');
    behavior_directory = strcat(session_directory, '/behavior/task/');
    
    %behavior_runs = struct2cell(dir(strcat(behavior_directory, '/l_full_*.dgz')));
    %runs = behavior_runs(1, :);
    
     bold_runs = struct2cell(dir(strcat(bold_directory, '0*')));
     runs = str2double(bold_runs(1, :));
%     
    for r = 1:length(runs)
        %this_run = str2double(runs{r}(end-5:end-4));
        %behavior_file = dir(fullfile(behavior_directory, strcat('*full*', sprintf('%02d', this_run), '.dgz')));
        behavior_file = dir(fullfile(behavior_directory, strcat('*full*', sprintf('%02d', runs(r)), '.dgz')));
        if length(behavior_file) > 0
            %finfo = dir(fullfile(bold_directory, sprintf('%03d', this_run), 'mot_demean*'));
            finfo = dir(fullfile(bold_directory, sprintf('%03d', runs(r)), 'mot_demean*'));
            mvt = load(fullfile(finfo.folder, finfo.name));
            run_indices = find(mvt(:, 1) ~= 0);
            
            censor_finfo = dir(fullfile(session_directory, 'preprocessed/motion_*censor.1D'));
            censor = load(fullfile(censor_finfo.folder, censor_finfo.name));
            run_censor = not(censor(run_indices));
            proportion = (sum(run_censor)/length(run_censor))*100;
            
            behavior_group = dg_read(fullfile(behavior_file.folder, behavior_file.name));
            v = string(behavior_group.variant);
            v = v(1);

            trials_to_analyze_rts = behavior_group.rts(3:end);
            trials_to_analyze_response = behavior_group.resp(3:end);
            trials_to_analyze_status = behavior_group.status;
            trials_to_analyze_side = behavior_group.side;
            
            percentcorrect = (round(mean(trials_to_analyze_status), 2))*100;
          
            sidebias = abs(round(mean(trials_to_analyze_status(trials_to_analyze_side == 0)), 2) - round(mean(trials_to_analyze_status(trials_to_analyze_side == 1)), 2));
            %nskiptrials = sum(trials_to_analyze_rts > 12000);
            nskiptrials = sum(trials_to_analyze_response == -1);
            %remaining_rts = trials_to_analyze_rts(trials_to_analyze_rts < 12000);
            remaining_rts = trials_to_analyze_rts(trials_to_analyze_response > -1);
            
            ntoolong = sum(remaining_rts >= ((prctile(remaining_rts, 75)) + (1.5 * iqr(remaining_rts))) & remaining_rts > 4000);
            
            if (percentcorrect >= 75) && (sidebias <= 0.25) && (nskiptrials <= 2) && (ntoolong <=3 ) && (proportion <= 15)
                use = 'Y';
            else
                use = 'N';
            end


            session = [session; sessions(i)];
            run = [run; runs(r)];
            variant = [variant; v];
            percent_correct = [percent_correct; percentcorrect];
            side_bias = [side_bias; sidebias];
            n_skip_trials = [n_skip_trials; nskiptrials];
            n_too_long = [n_too_long; ntoolong];
            percent_censor = [percent_censor; proportion];
            usable = [usable; use];
        end
    end
end
      
full_table = table(session, run, variant, percent_correct, side_bias, n_skip_trials, n_too_long, percent_censor, usable);
usable_table = full_table(full_table.usable == 'Y', :);
%usable_table = usable_table(usable_table.variant ~= "native", :);

perception_indices = find(usable_table.variant == "perception");
%usable_table([perception_indices(21:end)], :) = [];

control_indices = find(usable_table.variant == "counting");
%usable_table([control_indices(21:end)], :) = [];
        
        
    