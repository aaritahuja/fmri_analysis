load('SPM.mat')
contrasts = {'simulation*', 'perception*', 'counting*', 'native*', 'native_post_response*'};

%contrasts = {'native*'};
betas_good = [];
betas_bad = [];
rows_to_remove = [];

figure; hold on;

for i = 1:length(contrasts)
    all_correlations = [];
    locations = find(contains(SPM.xX.name, contrasts{i}));
    movement_columns = find(contains(SPM.xX.name, 'mvt'));
    
    %native_pre = find(contains(SPM.xX.name, 'native*'));
    %native_post = find(contains(SPM.xX.name, 'native_post_response*'));
    
    for l = 1:length(locations)
        run_correlations = [];
        current_column = locations(l);
        this_run_movement_columns = movement_columns(movement_columns > current_column);
        this_run_movement_columns = [this_run_movement_columns(1):this_run_movement_columns(1)+5];
        %movement_columns = native_post(l);
        current_rows = find(SPM.xX.X(:, current_column) ~= 0);
        
        betas_good = [betas_good; current_column];            
        %to_eliminate = 0;
        
        for m = 1:length(this_run_movement_columns)
            correlation = corrcoef(SPM.xX.X(min(current_rows):max(current_rows), current_column), SPM.xX.X(min(current_rows):max(current_rows), this_run_movement_columns(m)));
            all_correlations = [all_correlations; correlation(2)];
            run_correlations = [run_correlations; correlation(2)];
        end
        
        if max(abs(run_correlations)) >= 0.5
            disp(['Contrast ', num2str(contrasts{i}), ', run ', num2str(l), ' has a correlation of ', num2str(max(abs(run_correlations)))])
            betas_bad = [betas_bad; current_column];
            rows_to_remove = [rows_to_remove; l];
        end
        
    end
        
    
histogram(all_correlations, 'BinWidth', 0.1)
xlim([-1, 1])
ylim([0, 30])
mean(all_correlations)

end

legend('simulation', 'perception', 'control', 'native_pre', 'native_post')

%betas_good = betas_good(not(ismember(betas_good, betas_bad)));

% for i = 1:length(betas_good)
% current_file = niftiread(sprintf('beta_%04d', betas_good(i)));
% current_file = current_file(:);
% average_beta = mean(current_file, 'omitnan');
% good_averages = [good_averages; average_beta];
% end
% 
% for i = 1:length(betas_bad)
% current_file = niftiread(sprintf('beta_%04d', betas_bad(i)));
% current_file = current_file(:);
% average_beta = mean(current_file, 'omitnan');
% bad_averages = [bad_averages; average_beta];
% end
% 
% sessions = [5, 7, 10:13, 15:19];
% [final_table, usable_table] = run_table(MRI_directory, sessions);
% simulation_table = usable_table(usable_table.variant == "simulation", :);
% perception_table = usable_table(usable_table.variant == "perception", :);
% control_table = usable_table(usable_table.variant == "counting", :);
% native_table = usable_table(usable_table.variant == "native", :);
% simulation_table(18,:) = [];
% simulation_table(14,:) = [];
% simulation_table(10,:) = [];
% control_table(19,:) = [];
% control_table(15,:) = [];
% control_table(11:12,:) = [];
% control_table(8,:) = [];
% control_table(7,:) = [];
% control_table([1,3],:) = [];
% perception_table(14,:) = [];
% perception_table(9,:) = [];
% perception_table(8,:) = [];
% perception_table(7,:) = [];
% perception_table(5,:) = [];
% perception_table(1,:) = [];
% combined_table = [simulation_table_new; perception_table_new; control_table_new; native_table_new];
% combined_table = sortrows(combined_table,[1, 2]);