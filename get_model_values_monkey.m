%% Get model values (beta, standard error, and confidence interval) for each voxel for a given condition in an ROI
%MRI_directory is the name of where all MRI_data is stored
%subjid is the subject id
%ROI_file is the COMPLETE path to the ROI .nii saved from Marsbar (for example: 
%'/Users/Aaru/Documents/MRI_data/localizer_second_level/fixed_effects/9-Feb-2020/combined_motion.nii')
%modeltype is a string for the specifics of the model (for example:
%'task_without_duration')
%contrast is a string for the specific contrast to use (for example:
%'simulation_over_counting')
%Complete example:
%[simcountbeta, simcountSE, simcountCI] = get_model_values(MRI_directory, 12, 'task_without_duration', 'simulation_over_counting', ROI_file);

function [beta, t, SE, CI, all_coords, translated_coords, varargout] = get_model_values_monkey(model_directory, contrast, ROI_file)

%Get ROI voxel indices
ROI_epi = niftiread(ROI_file);
ROI = ROI_epi > 0;
[x, y, z] = ind2sub(size(ROI_epi), find(ROI == 1));
ROI_coordinates = [x y z];

% ROI = load(ROI_file);
%x = ROI.higher_order_motion_coordinates(:, 1);
%y = ROI.higher_order_motion_coordinates(:, 2);
%z = ROI.higher_order_motion_coordinates(:, 3);

beta = zeros(length(ROI_coordinates(:, 1)), 1);
SE = zeros(length(ROI_coordinates(:, 2)), 1);
CI = zeros(length(ROI_coordinates(:, 3)), 1);

%Get SPM file for model
cd(model_directory)
load('SPM.mat');

%Use ROI voxel indices to get model values for each index
for i = 1:length(beta)
    coords = [ROI_coordinates(i, 1); ROI_coordinates(i, 2); ROI_coordinates(i, 3)];
%   coords = [ROI_x(i); ROI_y(i); ROI_z(i)];
    trans_coords = SPM.xVol.M(1:3, :)*[coords; ones(1,size(coords,2))];
    %trans_coords = transpose(trans_coords);
    %[currbeta, currSE, currCI] = extractSPMData(SPM,trans_coords, contrast);
    [currbeta, currSE, currCI] = extractSPMData(SPM, coords, contrast);
%    fprintf('Voxel %d done \n', i)
    beta(i, 1) = currbeta;
    SE(i, 1) = currSE;
    CI(i, 1) = currCI;
    all_coords(i, :) = coords;
    translated_coords(i, :) = trans_coords;
end

t = rmmissing(beta./SE);

if nargout == 3
    varargout{1} = SE;
end

if nargout == 4
    varargout{2} = CI;
end