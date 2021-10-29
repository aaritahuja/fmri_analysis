function [all_RDMs, RDM_lists, average_RDM] = compute_similarity_monkey(structure, measurement_type, distance_type)

%figure
%hold on
% random_distance = 'n/a';
sessions = fieldnames(structure);
% real_distance.sim_per_rhos = zeros(length(subjects), 1);
% real_distance.sim_count_rhos = zeros(length(subjects), 1);
% real_distance.per_count_rhos = zeros(length(subjects), 1);

conditions = fieldnames(structure.(sessions{1}));

for s = 1:length(sessions)
    
    if strcmp(measurement_type, 't')
        %Get t values
        %native_pre = rmmissing(structure.(sessions{s}).(conditions{1}).t);
        %native_post = rmmissing(structure.(sessions{s}).(conditions{2}).t);
        sim = rmmissing(structure.(sessions{s}).(conditions{1}).t);
        per = rmmissing(structure.(sessions{s}).(conditions{2}).t);
        count = rmmissing(structure.(sessions{s}).(conditions{3}).t);
    elseif strcmp(measurement_type, 'beta')
        %Get beta values
        sim = rmmissing(structure.(sessions{s}).(conditions{1}).beta);
        per = rmmissing(structure.(sessions{s}).(conditions{2}).beta);
        count = rmmissing(structure.(sessions{s}).(conditions{3}).beta);
    else
        disp('Measurement type not recognized. t or beta are acceptable types')
    end
    
    %Could be returned as function outputs if desired
    mean_sim(s) = mean(sim);
    mean_per(s) = mean(per);
    mean_count(s) = mean(count);
    
    %combined_matrix = [transpose(native_pre); transpose(native_post); transpose(sim); transpose(per); transpose(count)];
    combined_matrix = [transpose(sim); transpose(per); transpose(count)];
    
    if strcmp(distance_type, 'euclidean')
        %Get euclidean distance
        distance = pdist(combined_matrix, 'euclidean');
        distance = distance/length(sim);
    elseif strcmp(distance_type, 'correlation')
        %Get correlation values
        distance = pdist(combined_matrix, 'spearman');
    end
     
     RDM_lists(:, :, s) = distance;
     RDM_matrix(:, :, s) = squareform(distance);
end

all_RDMs = RDM_matrix;
average_RDM = mean(all_RDMs, 3);
figure
heat = heatmap(average_RDM);
labels = {'Simulation', 'Perception', 'Control'};
heat.XDisplayLabels = labels;
heat.YDisplayLabels = labels;
heat.NodeChildren(3).YDir='normal'; 

    