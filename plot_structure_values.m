function [beta_matrix, beta_mean] = plot_structure_values(structure, model_type, contrasts, ROI_file, normalize_scale)


ROI_file = '/Users/Aaru/Documents/MRI_data/localizer_second_level/motion_over_static/23-Mar-2020/higher_order_motion.mat';
ROI = load(ROI_file);
x = ROI.higher_order_motion_coordinates(:, 1);
y = ROI.higher_order_motion_coordinates(:, 2);
z = ROI.higher_order_motion_coordinates(:, 3);

%Set some basic info, subjects, ROI etc.
% subjects = fieldnames(structure);
 subjects = {'s3'};
% ROI_epi = readnifti(ROI_file);
% ROI = ROI_epi > 0;
% [x, y, z] = ind2sub(size(ROI_epi), find(ROI == 1));

meta_beta_matrix = zeros(length(z), length(contrasts));
%meta_CI_matrix = zeros(length(z), length(contrasts));

for c = 1:length(contrasts)    
    beta_matrix = zeros(length(structure.s3.(model_type).counting.beta), length(subjects));
    %CI_matrix = zeros(length(structure.s2.(model_type).counting.beta), length(fieldnames(structure)));
    for i = 1:length(subjects)
        beta_matrix(:, i) = structure.(subjects{i}).(model_type).(contrasts{c}).beta;
        %CI_matrix(:, i) = structure.(subjects{i}).(model_type).(contrasts{c}).CI;
    end
    beta_mean = nanmean(beta_matrix, 2);
    %CI_mean = nanmean(CI_matrix, 2);
    
    meta_beta_matrix(:, c) = beta_mean;
    %meta_CI_matrix(:, c) = CI_mean;
end

%Get values to later normalize all data to same color scale
colormax = max(max(meta_beta_matrix));
colormin = min(min(meta_beta_matrix));

%Normalize all CI values to same size scale
%meta_CI_matrix = ((max(max(meta_CI_matrix)) + 1) - meta_CI_matrix);
%meta_CI_matrix = (30/max(max(meta_CI_matrix)))*meta_CI_matrix;

%draw
%close(gcf)
draw = figure;
for c = 1:length(contrasts)
    subplot(1, length(contrasts), c)
    if normalize_scale == 1
        scatter3(x, y, z, meta_CI_matrix(:, c), meta_beta_matrix(:, c), 'filled')
    else
        scatter3(x, y, z, 20, meta_beta_matrix(:, c), 'filled')
    end
    view(54, 50)
    title(contrasts{c})
    %colormap heat
    caxis([colormin colormax])
    set(gca,'visible','off')    
    axis vis3d
end
figure(draw)
set(gcf, 'Position', [500 1000 1600 500])