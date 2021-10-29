simulation_all = [];
perception_all = [];
control_all = [];

for r = 1:height(usable_table)
    id_number = '%03d';
    session_number = sprintf(id_number, usable_table.session(r));
    run_number = sprintf(id_number, usable_table.run(r));
    current_directory = fullfile(MRI_directory, strcat('Glenn/sub-glenn', session_number), 'ses-01/bold', run_number);
    cd(current_directory)
    
    files = dir('s3wrlf*');
    all_files = {files.name};
    
    for i = 1:length(all_files)
        current_file = niftiread(all_files{i});
        current_file = current_file(:);
        current_file = current_file(mask);
        if usable_table.variant(r) == "simulation"
            simulation_all = [simulation_all, current_file];
        elseif usable_table.variant(r) == "perception"
            perception_all = [perception_all, current_file];
        elseif usable_table.variant(r) == "counting"
            control_all = [control_all, current_file];
        end
        %disp(i)
    end
    
    disp(strcat('Session ', session_number, ', ', 'run ', run_number, ' done'))
    disp(strcat(num2str(height(usable_table) - r), ' runs left'))
end

simulation_averaged = mean(simulation_all, 2);
perception_averaged = mean(perception_all, 2);
control_averaged = mean(control_all, 2);
histogram(simulation_averaged); hold on
histogram(perception_averaged);
histogram(control_averaged);
legend('simulation', 'perception', 'control')