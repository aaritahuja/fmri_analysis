function [average_simulation_densities, average_perception_densities, average_counting_densities] = check_saccade_densities(MRI_directory, subjids, task_runs)

for i = 1:length(subjids)
    subjid = subjids(i);
    task_run = cell2mat(task_runs(i));
    n_runs = length(task_run);
    subject_directory = fullfile(MRI_directory, num2str(subjid));
    behavior_directory = fullfile(subject_directory, 'behavior/task');
    simulation_densities = [];
    counting_densities = [];
    perception_densities = [];
    for r = 1:n_runs
        behavior_file = dir(fullfile(behavior_directory, strcat('*0', num2str(task_run(r)), '.dgz')));
        behavior_group = dg_read(fullfile(behavior_directory, behavior_file.name));
        variant_column = cellstr(behavior_group.variant);
        variant = variant_column{1};
        for t = 1:length(behavior_group.fixon)
            n_valid_sacs = sum(behavior_group.stimon(t) < behavior_group.sactimes{t, 1} &  behavior_group.sactimes{t, 1} < behavior_group.response(t));
            sac_density = n_valid_sacs/(behavior_group.rts(t)/1000);
            run_density(t, 1) = sac_density;
        end
        mean_run_density = mean(run_density);
        if strcmp(variant, 'simulation')
            simulation_densities(end+1) = mean_run_density;
        elseif strcmp(variant, 'counting')
            counting_densities(end+1) = mean_run_density;
        elseif strcmp(variant, 'perception')
            perception_densities(end+1) = mean_run_density;
        end
    end
    average_simulation_densities(i, 1) = mean(simulation_densities);
    average_perception_densities(i, 1) = mean(perception_densities);
    average_counting_densities(i, 1) = mean(counting_densities);
end