function prep_timing_file_afni(MRI_directory, output_directory, usable_table, censor_bad_trials)

simulation = [];
perception = [];
control = [];
master_skip = [];
master_outlier = [];
master_initial = [];
combined_censor = [];
cumulative_files = 0
cumulative_censor = 0

for r = 1:length(usable_table.session)
    sess = sprintf('%03d', usable_table.session(r));
    session_directory = fullfile(MRI_directory, strcat('Glenn/sub-glenn', sess), 'ses-01');
    behavior_directory = strcat(session_directory, '/behavior/task/');
    mkdir(fullfile(behavior_directory, 'usable_new'));
    behavior_file = dir(fullfile(behavior_directory, strcat('*full*', sprintf('%02d', usable_table.run(r)), '.dgz')));
    behavior_group = dg_read(fullfile(behavior_file.folder, behavior_file.name));
    copyfile(fullfile(behavior_directory, behavior_file.name), fullfile(behavior_directory, 'usable_new'));
    
    if censor_bad_trials == 1
%%
        trials = [];
        too_long_scan_indices = [];
        skip_trial_scan_indices = [];
        
        %Get trial onsets and durations
        for t = 1:length(behavior_group.pre_resp_durations)
            trials = [trials, strcat(num2str(behavior_group.pre_resp_onsets(t)), ":", num2str(behavior_group.pre_resp_durations(t)))];
        end
        
        placeholder = [repmat("*", 1, length(trials))];
        v = string(behavior_group.variant);
        v = v(1);
        
        %Organize onsets and durations into tables based on the variant
        if v == "simulation"
            simulation = [simulation; trials];
        else
            simulation = [simulation; placeholder];
        end
        
        if v == "perception"
            perception = [perception; trials];
        else
            perception = [perception; placeholder];
        end
        
        if v == "counting"
            control = [control; trials];
        else
            control = [control; placeholder];
        end
        
        tr = 1.8; 
        
        %Get scan indices for scans spanning the start of the first trial to the start of the third trial
        first_two_trials = behavior_group.pre_resp_onsets([1, 3]);
        first_scan_times = first_two_trials/tr;
        first_scan_indices = transpose([round(first_scan_times(1)):floor(first_scan_times(2) - 1)]);
        
        %Get scan indices for trials with outlier reaction times
        too_long = find(behavior_group.rts(3:end) >= ((prctile(behavior_group.rts(3:end), 75)) + (1.5 * iqr(behavior_group.rts(3:end)))) & behavior_group.rts(3:end) > 4000);
        too_long = too_long + 2;
        for tl = 1:length(too_long)
            too_long_trial_start = behavior_group.pre_resp_onsets(too_long(tl));
            %too_long_trial_end = behavior_group.pre_resp_onsets(too_long(tl)) + (behavior_group.stimoff(too_long(tl))/1000);
            too_long_trial_end = behavior_group.pre_resp_onsets(too_long(tl)) + (behavior_group.stimoff(too_long(tl))/1000) + (behavior_group.feedback_iti(too_long(tl))/1000);
            %too_long_trial_end = behavior_group.pre_resp_onsets(too_long(tl) + 1);
            too_long_scans = transpose([round((too_long_trial_start)/tr):floor((too_long_trial_end)/1.8)]);
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
        
        response_times = behavior_group.response_onsets;
        response_scans = round(behavior_group.response_onsets/tr);
        
        %Get number of scans in this run by looking in the bold directory (this
        %is temporary, come up with a better system)
        bold_directory = fullfile(session_directory, 'bold', sprintf('%03d', usable_table.run(r)));
        bold_finfo = dir(fullfile(bold_directory, '*.nii'));
        n_files = length(bold_finfo);
        cumulative_files = cumulative_files + n_files
        
        %Get the movement regressor file to derive the indices that correspond
        %to this run
        preprocessed_directory = fullfile(session_directory, 'preprocessed');
        finfo = dir(fullfile(preprocessed_directory, strcat('mot_demean.r', sprintf('%02d', usable_table.run(r)), '.1D')));
        mvt = load(fullfile(finfo.folder, finfo.name));
        run_indices = find(mvt(:, 1) ~= 0);
        
        %Use these run indices to select the appropriate volumes to censor
        censor_finfo = dir(fullfile(preprocessed_directory, 'motion_*censor.1D'));
        censor = load(fullfile(censor_finfo.folder, censor_finfo.name));
        run_censor = censor(run_indices);
        
        %Manually add in censoring for the first to trials, skip trials, and
        %outlier rt trials. Then concat this array to cumulative array
        run_censor([skip_trial_scan_indices; too_long_scan_indices; first_scan_indices; response_scans]) = 0;
        cumulative_censor = cumulative_censor + sum(not(run_censor))
        combined_censor = [combined_censor; run_censor];
        
    elseif censor_bad_trials == 0
%%
        skip = [];
        outlier = [];
        initial = [];
        trials = [];
        
        first_two = find(behavior_group.pre_resp_onsets(1:2));
        
        too_long = find(behavior_group.rts(3:end) >= ((prctile(behavior_group.rts(3:end), 75)) + (1.5 * iqr(behavior_group.rts(3:end)))) & behavior_group.rts(3:end) > 4000);
        too_long = too_long + 2;
        
        skip_trials = find(behavior_group.rts(3:end) > 12000);
        skip_trials = skip_trials + 2;
                 
        remaining = find(behavior_group.pre_resp_onsets);
        remaining_selector = ones(size(remaining));
        remaining_selector([first_two; skip_trials; too_long]) = 0;
        remaining = remaining(logical(remaining_selector));
        
                
        for t = 1:length(behavior_group.pre_resp_durations)
            if any(t == too_long)
                outlier = [outlier, strcat(num2str(behavior_group.pre_resp_onsets(t)), ":", num2str(behavior_group.pre_resp_durations(t)))];
            else
                outlier = [outlier, "*" ];
            end
            
            if any(t == skip_trials)
                skip = [skip, strcat(num2str(behavior_group.pre_resp_onsets(t)), ":", num2str(behavior_group.pre_resp_durations(t)))];
            else
                skip = [skip, "*"];
            end
            
            if any(t == first_two)
                initial = [initial, strcat(num2str(behavior_group.pre_resp_onsets(t)), ":", num2str(behavior_group.pre_resp_durations(t)))];
            else
                initial = [initial, "*"];
            end
            
            if any(t == remaining)
                trials = [trials, strcat(num2str(behavior_group.pre_resp_onsets(t)), ":", num2str(behavior_group.pre_resp_durations(t)))];
            else
                trials = [trials, "*"];
            end
        end
        
        placeholder = [repmat("*", 1, length(trials))];
        v = string(behavior_group.variant);
        v = v(1);
        
        if v == "simulation"
            simulation = [simulation; trials];
        else
            simulation = [simulation; placeholder];
        end
        
        if v == "perception"
            perception = [perception; trials];
        else
            perception = [perception; placeholder];
        end
        
        if v == "counting"
            control = [control; trials];
        else
            control = [control; placeholder];
        end
        
        master_skip = [master_skip; skip];
        master_outlier = [master_outlier; outlier];
        master_initial = [master_initial; initial];
        
        preprocessed_directory = fullfile(session_directory, 'preprocessed');
        finfo = dir(fullfile(preprocessed_directory, strcat('mot_demean.r', sprintf('%02d', usable_table.run(r)), '.1D')));
        mvt = load(fullfile(finfo.folder, finfo.name));
        run_indices = find(mvt(:, 1) ~= 0);
        
        %Use these run indices to select the appropriate volumes to censor
        censor_finfo = dir(fullfile(preprocessed_directory, 'motion_*censor.1D'));
        censor = load(fullfile(censor_finfo.folder, censor_finfo.name));
        run_censor = censor(run_indices);
        cumulative_censor = cumulative_censor + sum(not(run_censor))
        combined_censor = [combined_censor; run_censor];
        
        %Get number of scans in this run by looking in the bold directory (this
        %is temporary, come up with a better system)
        bold_directory = fullfile(session_directory, 'bold', sprintf('%03d', usable_table.run(r)));
        bold_finfo = dir(fullfile(bold_directory, '*.nii'));
        n_files = length(bold_finfo);
        cumulative_files = cumulative_files + n_files

    end
end
    
simulation_table = table(simulation);
perception_table = table(perception);
control_table = table(control);
sum(not(combined_censor))

censor_table = table(combined_censor);

if censor_bad_trials == 0
    initial_table = table(master_initial);
    skip_table = table(master_skip);
    outlier_table = table(master_outlier);
end

out = fullfile(MRI_directory, 'Glenn/task_models/afni', output_directory);

if ~exist(out, 'dir')
    mkdir(out)
end

cd(out)
writetable(simulation_table, 'simulation_timings.txt', 'WriteVariableNames', 0, 'Delimiter', ' ');
writetable(perception_table, 'perception_timings.txt', 'WriteVariableNames', 0, 'Delimiter', ' ');
writetable(control_table, 'control_timings.txt', 'WriteVariableNames', 0, 'Delimiter', ' ');
writetable(censor_table, 'censor_indices.txt', 'WriteVariableNames', 0, 'Delimiter', ' ');

if censor_bad_trials == 0
    writetable(initial_table, 'first_trials.txt', 'WriteVariableNames', 0, 'Delimiter', ' ');
    writetable(skip_table, 'skip_trials.txt', 'WriteVariableNames', 0, 'Delimiter', ' ');
    writetable(outlier_table, 'outlier_trials.txt', 'WriteVariableNames', 0, 'Delimiter', ' ');
end
end