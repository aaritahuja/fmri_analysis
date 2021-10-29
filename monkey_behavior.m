%Example: master_list = glenn_behavior('/Users/Aaru/Documents/MRI_monkey', [5:7, 9:14])

function master_list = glenn_behavior(MRI_directory, behavior_table)


performance = [];
uncertainty = [];
rts = [];
intersection = [];
variant = [];
plank_angle_prediction = [];
plank_clustering_prediction = [];
valid_first_bounce = [];
file_list = [];
session_list = [];
run_list = [];
remove = [];

shuffled_intersection = [];

for i = 1:height(behavior_table)
    session = sprintf('%03d', behavior_table.session(i));
    behavior_directory = strcat(MRI_directory, '/Glenn/sub-glenn', session, '/ses-01/behavior/task/');
    %behav_files = struct2cell(dir(strcat(behavior_directory, '*full*')));
    %behav_files = behav_files(1, :);
    run = behavior_table.run(i);
    behavior_file = dir(fullfile(behavior_directory, strcat('*full*', sprintf('%02d', run), '.dgz')));
    behavior_group = dg_read(fullfile(behavior_file.folder, behavior_file.name));
    
    %trials_to_analyze_rts = behavior_group.rts(3:end);
    remaining_rts = behavior_group.rts(3:end);
    remaining_resp = behavior_group.resp(3:end);
    
    %skip_trials = trials_to_analyze_rts > 12000;
    skip_trials = remaining_resp == -1;
    skip_idx = find(skip_trials);
    
    rts_to_check = remaining_rts(remaining_resp > -1);
    
    long_trials = rts_to_check >= ((prctile(rts_to_check, 75)) + (1.5 * iqr(rts_to_check))) & rts_to_check > 4000;
    if length(skip_idx) > 0
        for s = 1:length(skip_idx)
            current_idx = skip_idx(s);
            long_trials = [long_trials(1:current_idx-1); 0; long_trials(current_idx:end)];
        end
    end
    
    trials_to_remove = long_trials + skip_trials;
    %trials_to_remove = skip_trials;
    remove = [remove; 1; 1; trials_to_remove];
    
    performance = [performance; behavior_group.status];
    uncertainty = [uncertainty; behavior_group.discrete_uncertainty];
    rts = [rts; behavior_group.rts];
    %intersection = [intersection; behavior_group.intersection];
    %shuffled_intersection = [shuffled_intersection; repelem(mean(behavior_group.shuffled_intersections{1, 1}), length(behavior_group.intersection), 1)];
    variant = [variant; string(behavior_group.variant)];
    plank_angle_prediction = [plank_angle_prediction; behavior_group.plank_angle_prediction];
    plank_clustering_prediction = [plank_clustering_prediction; behavior_group.plank_clustering_prediction];
    valid_first_bounce = [valid_first_bounce; behavior_group.trajectory_prediction_valid];
    file_list = [file_list; string(repmat(behavior_file.name, length(behavior_group.variant), 1))];
    session_list = [session_list; repmat(session, length(behavior_group.variant), 1)];
    run_list = [run_list; repmat(run, length(behavior_group.variant), 1)];
    

end

%master_list = table(file_list, session_list, variant, uncertainty, rts, performance, intersection, shuffled_intersection, plank_angle_prediction, plank_clustering_prediction, valid_first_bounce, remove);
master_list = table(file_list, session_list, run_list, variant, uncertainty, rts, performance, plank_angle_prediction, plank_clustering_prediction, valid_first_bounce, remove);
master_list = master_list((remove== 0), :);

%mean_intersections = varfun(@mean, master_list, 'GroupingVariables', {'file_list'}, 'InputVariables', {'intersection', 'variant'});
%intersection_comparison = table(mean_intersections.variant, shuffled_intersection)

test = varfun(@mean, master_list, 'GroupingVariables', {'variant', 'uncertainty'}, 'InputVariables', {'performance', 'rts'});
