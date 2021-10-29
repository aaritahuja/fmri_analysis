function master_list = human_behavior(MRI_directory, subjids)
subject = [];
variant = [];
balldist = [];
nhit = [];
uncertainty = [];
performance = [];
rts = [];


for s = 1:length(subjids)
    subjid = subjids(s);
    subject_directory = fullfile(MRI_directory, num2str(subjid));
    behavior_directory = fullfile(subject_directory, 'behavior/task');
    
    conditions = {'simulation' 'perception' 'counting'};
    for c = 1:length(conditions)
        behavior_files = struct2cell(dir(fullfile(behavior_directory, strcat('*', conditions{c}, '*.dgz'))));
        behavior_files = behavior_files(1, :);
        for i = 1:length(behavior_files)
            behavior_file = behavior_files{i};
            behavior_group = dg_read(fullfile(behavior_directory, behavior_file));
            %trials_to_remove = behavior_group.rts >= ((prctile(behavior_group.rts, 75)) + (1.5 * iqr(behavior_group.rts))) & behavior_group.rts > 3000;
            trials_to_remove = repmat(0, length(behavior_group.rts), 1);
            subject = [subject; repmat(subjid, (length(behavior_group.rts)-sum(trials_to_remove)), 1)];
            
            this_run_variant = string(behavior_group.variant);
            variant = [variant; this_run_variant(not(trials_to_remove))];
            
            balldist = [balldist; behavior_group.balldist(not(trials_to_remove))];
            nhit = [nhit; behavior_group.nhit(not(trials_to_remove))];
            uncertainty = [uncertainty; behavior_group.discrete_uncertainty(not(trials_to_remove))];
            performance = [performance; behavior_group.status(not(trials_to_remove))];
            rts = [rts; behavior_group.rts(not(trials_to_remove))];
            
            %variant = [variant; string(behavior_group.variant)];
            %variant = [variant; string(behavior_group.variant)];
%% Eye movement stuff
%            for t = 1:length(behavior_group.stimon)
%                 n_valid_sacs = sum(behavior_group.stimon(t) < behavior_group.sactimes{t, 1} &  behavior_group.sactimes{t, 1} < behavior_group.response(t));
%                 sac_density = n_valid_sacs/(behavior_group.rts(t)/1000);
%                 sac_density_counter(t, 1) = sac_density;
%                 sac_counter(t, 1) = n_valid_sacs;
%            end
%             too_long = not(behavior_group.rts >= ((prctile(behavior_group.rts, 50)) + (1.5 * iqr(behavior_group.rts))));
%             sac_counter = sac_counter(too_long);
%             sac_density_counter = sac_density_counter(too_long);
%             rts = rts(too_long);
%             average_sacs = mean(sac_counter);
%             average_sac_density = mean(sac_density_counter);
%             average_rt = mean(rts);
%             metagroup(i, 1) = average_sacs;
%             metagroup(i, 2) = average_sac_density;
%             metagroup(i, 3) = average_rt;
        end
%         mean_low_uncertainty_rt = mean(rts(uncertainty == 0));
%         mean_high_uncertainty_rt = mean(rts(uncertainty == 1));
%         mean_low_uncertainty_performance = mean(performance(uncertainty == 0));
%         mean_high_uncertainty_performance = mean(performance(uncertainty == 0));
%         mean_short_rt = mean(rts(nhit == 2));
%         mean_long_rt = mean(rts(nhit == 4));
%         mean_short_performance = mean(performance(nhit == 2));
%% More eye movement stuff        
%         average_sacs_this_condition = mean(metagroup(:, 1));
%         average_sac_density_this_condition = mean(metagroup(:, 2));
%         average_rt_this_condition = mean(metagroup(:, 3));
%         if strcmp(conditions{c}, 'simulation')
%             simulation_sac_density(s, 1) = average_sac_density_this_condition;
%             simulation_sacs(s, 1) = average_sacs_this_condition;
%             simulation_rts(s, 1) = average_rt_this_condition;
%         elseif strcmp(conditions{c}, 'perception')
%             perception_sac_density(s, 1) = average_sac_density_this_condition;
%             perception_sacs(s, 1) = average_sacs_this_condition;
%             perception_rts(s, 1) = average_rt_this_condition;
%         elseif strcmp(conditions{c}, 'counting')
%             counting_sac_density(s, 1) = average_sac_density_this_condition;
%             counting_sacs(s, 1) = average_sacs_this_condition;
%             counting_rts(s, 1) = average_rt_this_condition;
%         end
    end

end
rts = rts/1000;
master_list = table(subject, variant, balldist, nhit, uncertainty, performance, rts);
%test = varfun(@mean, master_list_2, 'GroupingVariables',{'balldist'});