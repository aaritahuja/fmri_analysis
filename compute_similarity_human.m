function [all_RDMs, RDM_lists, average_RDM] = compute_similarity (structure, model_type, measurement_type, distance_type, varargin)

%figure
%hold on
% random_distance = 'n/a';
subjects = fieldnames(structure);
% real_distance.sim_per_rhos = zeros(length(subjects), 1);
% real_distance.sim_count_rhos = zeros(length(subjects), 1);
% real_distance.per_count_rhos = zeros(length(subjects), 1);

conditions = struct2cell(structure.(subjects{1}));
conditions = conditions{1, 1};
conditions = fieldnames(conditions);

for s = 1:length(subjects)
    
    if strcmp(measurement_type, 't')
        %Get t values
        sim = rmmissing(structure.(subjects{s}).(model_type).(conditions{1}).t);
        per = rmmissing(structure.(subjects{s}).(model_type).(conditions{2}).t);
        count = rmmissing(structure.(subjects{s}).(model_type).(conditions{3}).t);
    elseif strcmp(measurement_type, 'beta')
        %Get beta values
        sim = rmmissing(structure.(subjects{s}).(model_type).(conditions{1}).beta);
        per = rmmissing(structure.(subjects{s}).(model_type).(conditions{2}).beta);
        count = rmmissing(structure.(subjects{s}).(model_type).(conditions{3}).beta);
    else
        disp('Measurement type not recognized. t or beta are acceptable types')
    end
    
    %Could be returned as function outputs if desired
    mean_sim(s) = mean(sim);
    mean_per(s) = mean(per);
    mean_count(s) = mean(count);
    
    combined_matrix = [transpose(sim); transpose(per); transpose(count)];
    
    if strcmp(distance_type, 'euclidean')
        %Get euclidean distance
        distance = pdist(combined_matrix, 'euclidean');
        distance = distance/length(sim);
    elseif strcmp(distance_type, 'correlation')
        %Get correlation values
        distance = pdist(combined_matrix, 'spearman');
        distance = 1-distance;
    elseif strcmp(distance_type, 'mahalanobis');
        subject = reshape(subjects{s}, 1, [])';
        if length(subject) == 2
            subject = subject(2);
        elseif length(subject) == 3
            subject = strcat(subject(2), subject(3));
        end
        
        cd(fullfile('/Users/Aaru/Documents/MRI_data', subject, 'models', model_type ));
        %          %Get residual files
        residual_files = dir('Res_*.nii');
        residual_files = struct2cell(residual_files);
        residual_files(2:6, :) = [];
        %%
        %          %Extract residuals from ROI
        %          if length(varargin) < 1
        %              error('ROI_file not specified - not sure where to grab residuals from')
        %          else
        %              ROI_file = varargin{1};
        %          end
        %          %temp
        %          ROI_file = fullfile(ROI_file, subject, 'models', model_type, '/V1_resliced');
        %
        %          ROI_epi = readnifti(ROI_file);
        %          ROI = ROI_epi > 0;
        %          [x, y, z] = ind2sub(size(ROI_epi), find(ROI == 1));
        %%
%         indices_to_keep_a = find(isnan(structure.(subjects{s}).(model_type).(conditions{1}).beta) == 0);
%         indices_to_keep_b = find(structure.(subjects{s}).(model_type).(conditions{1}).beta ~= 0);
%         indices_to_keep = intersect(indices_to_keep_a, indices_to_keep_b);
%         sim = rmmissing(structure.(subjects{s}).(model_type).(conditions{1}).beta(indices_to_keep));
%         per = rmmissing(structure.(subjects{s}).(model_type).(conditions{2}).beta(indices_to_keep));
%         count = rmmissing(structure.(subjects{s}).(model_type).(conditions{3}).beta(indices_to_keep));
%         
%         x = structure.(subjects{s}).(model_type).(conditions{1}).coords(indices_to_keep, 1);
%         y = structure.(subjects{s}).(model_type).(conditions{1}).coords(indices_to_keep, 2);
%         z = structure.(subjects{s}).(model_type).(conditions{1}).coords(indices_to_keep, 3);

         x = structure.(subjects{s}).(model_type).(conditions{1}).coords(:, 1);
         y = structure.(subjects{s}).(model_type).(conditions{1}).coords(:, 2);
         z = structure.(subjects{s}).(model_type).(conditions{1}).coords(:, 3);
         
         %Make T timepoints X P voxels noise matrix
         noise_matrix = [];
         for t = 1:length(residual_files)
             residual = readnifti(residual_files{t});
             for i = 1:length(x)
                 noise_matrix(t, i) = [residual(x(i), y(i), z(i))];
             end
             %disp(strcat(residual_files{t}, ' done'))
         end
         
         noise_matrix = rmmissing(noise_matrix, 2);
         
         %Make P X P variance-covariance matrix
         covariance_matrix = (transpose(noise_matrix) * noise_matrix) / length(noise_matrix(:, 1));
         %nearest_covariance = nearestSPD(covariance_matrix);
         
         %Shrink the covariance matrix to make it invertible
         h = 0.4;
         shrunk = h*diag(diag(covariance_matrix)) + (1 - h)*covariance_matrix;
         %shrunk = covariance_matrix;

         pre_whitened_sim = transpose(sim) * (shrunk^-0.5);
         %sim_SE = rmmissing(transpose(structure.(subjects{s}).(model_type).(conditions{1}).SE ));
         %pre_whitened_sim = pre_whitened_sim./sim_SE;

         pre_whitened_per = transpose(per) * (shrunk^-0.5);
         %per_SE = rmmissing(transpose(structure.(subjects{s}).(model_type).(conditions{2}).SE ));
         %pre_whitened_per = pre_whitened_per./per_SE;

         pre_whitened_count = transpose(count) * (shrunk^-0.5);
         %count_SE = rmmissing(transpose(structure.(subjects{s}).(model_type).(conditions{3}).SE ));
         %pre_whitened_count = pre_whitened_count./count_SE;
         
         %Error is j
         
         %Inverse variance-covariance matrix (with test to make sure that
         %it is correct)
         %inverse_covariance_matrix = inv(shrunk);
         %test = inverse_covariance_matrix * shrunk;
         %inverse_covariance_matrix = inv(nearest_covariance);
         %test = inverse_covariance_matrix * nearest_covariance;
         %test = round(test);
         %if max(max(test)) ~= 1 || min(min(test)) ~= 0
         %    warning(strcat('Inverse covariance matrix for subject ', subject, ' might not be correct'))
         %end
         
         %calculate Mahalanobis distance as the square root of A-B * inverse covariance matrix * transposed A-B
         %sim_per_distance = sqrt(transpose((sim - per)) * inverse_covariance_matrix * (sim-per));
         %sim_count_distance = sqrt(transpose((sim - count)) * inverse_covariance_matrix * (sim-count));
         %per_count_distance = sqrt(transpose((per - count)) * inverse_covariance_matrix * (per-count));
         
         sim_per_distance = sqrt((pre_whitened_sim - pre_whitened_per) * transpose(pre_whitened_sim - pre_whitened_per));
         %se = std(pre_whitened_sim - pre_whitened_per)/sqrt(length(pre_whitened_sim));
         %sim_per_distance = sim_per_distance/se;
         
         sim_count_distance = sqrt((pre_whitened_sim - pre_whitened_count) * transpose(pre_whitened_sim - pre_whitened_count));
         %se = std(pre_whitened_sim - pre_whitened_count)/sqrt(length(pre_whitened_sim));
         %sim_count_distance = sim_count_distance/se;
         
         per_count_distance = sqrt((pre_whitened_per - pre_whitened_count) * transpose(pre_whitened_per - pre_whitened_count));
         %se = std(pre_whitened_per - pre_whitened_count)/sqrt(length(pre_whitened_per));
         %per_count_distance = per_count_distance/se;
         
         distance = [sim_per_distance, sim_count_distance, per_count_distance];
         
         disp(strcat('Subject ', subject, ' done'));
%%         
%          covariance = cov(combined_matrix);
%          SPD = nearestSPD(covariance);
%          inverse_SPD = inv(SPD);
%          sim_per_distance = (combined_matrix(1, :) - combined_matrix(2, :)) * inverse_SPD * (transpose(combined_matrix(1, :) - combined_matrix(2, :)));
%          sim_count_distance = (combined_matrix(1, :) - combined_matrix(3, :)) * inverse_SPD * (transpose(combined_matrix(1, :) - combined_matrix(3, :)));
%          per_count_distance = (combined_matrix(2, :) - combined_matrix(3, :)) * inverse_SPD * (transpose(combined_matrix(2, :) - combined_matrix(3, :)));
%          distance = [sim_per_distance, sim_count_distance, per_count_distance];
     end
     
     RDM_lists(:, :, s) = distance;
     RDM_matrix(:, :, s) = squareform(distance);

%%     
    %Plot
%     plot([sim_per_distance, sim_count_distance, per_count_distance], '-o')
%     xlim([0, 4])
%     xlab = cellstr(["", "Simulation-Perception", "Simulation-Counting", "Perception-Counting", ""]);
%     set(gca,'XTick',0:4);
%     set(gca,'XtickLabel',xlab(1:5));
    
%     %Get chance correlation
%     random_distance.random_sim_per_rhos.(strcat('s', num2str(s))) = zeros(100, 1);
%     random_distance.random_sim_count_rhos.(strcat('s', num2str(s))) = zeros(100, 1);
%     random_distance.random_per_count_rhos.(strcat('s', num2str(s))) = zeros(100, 1);
%     
%     for n = 1:100      
%         %Shuffle_lists
%         random_sim = rmmissing(structure.(subjects{s}).task_with_duration.simulation.t(randperm(size(structure.(subjects{s}).task_with_duration.simulation.t, 1)), :));
%         random_per = rmmissing(structure.(subjects{s}).task_with_duration.perception.t(randperm(size(structure.(subjects{s}).task_with_duration.perception.t, 1)), :));
%         random_count = rmmissing(structure.(subjects{s}).task_with_duration.counting.t(randperm(size(structure.(subjects{s}).task_with_duration.counting.t, 1)), :));
%         
%         %Get correlation values
%         [rho_sim_per_rand, pval_sim_per_rand] = corr(random_sim, random_per);
%         [rho_sim_count_rand, pval_sim_count_rand] = corr(random_sim, random_count);
%         [rho_per_count_rand, pval_per_count_rand] = corr(random_per, random_count);
%         
%         %Get euclidean distance values
% %         rho_sim_per_rand = mean(abs(random_sim - random_per));
% %         rho_sim_count_rand = mean(abs(random_sim - random_count));
% %         rho_per_count_rand = mean(abs(random_per - random_count));
%         
%         %Save values
%         random_distance.random_sim_per_rhos.(strcat('s', num2str(s)))(n) = rho_sim_per_rand;
%         random_distance.random_sim_count_rhos.(strcat('s', num2str(s)))(n) = rho_sim_count_rand;
%         random_distance.random_per_count_rhos.(strcat('s', num2str(s)))(n) = rho_per_count_rand;
%     end  
end

all_RDMs = RDM_matrix;
average_RDM = mean(all_RDMs, 3);
figure
heatmap(average_RDM);

    