function [first_map, second_map] = compare_maps(model_directory, map_1_label, map_1_path, map_2_label, map_2_path)
    
    
    mask = niftiread(fullfile(model_directory, 'resampled_mask.nii'));
    %mask = niftiread(fullfile(model_directory, 'motion_over_flicker_FWE.nii'));
    map_1_path = fullfile(model_directory, map_1_path);
    map_2_path = fullfile(model_directory, map_2_path);
    
    mask = mask(:);
    mask = (mask ~= 0);

    map_1 = niftiread(map_1_path);
    map_1 = map_1(:);
    map_1 = map_1(mask);
    
    map_2 = niftiread(map_2_path);
    map_2 = map_2(:);
    map_2 = map_2(mask);
    
    first_map = map_1;
    second_map = map_2;
    
    map_1_nonzero = map_1 ~= 0;
    map_2_nonzero = map_2 ~=0;
    combined = (map_1_nonzero + map_2_nonzero) == 2;
    map_1 = map_1(combined);
    map_2 = map_2(combined);
    
    mean(map_1, 'omitnan');
    mean(map_2, 'omitnan');
    
    figure; scatter(map_1, map_2, 1, 'red')
    xlabel(map_1_label)
    ylabel(map_2_label)
    title(strcat(map_1_label, 'vs', map_2_label))
    
    axis_min = floor(min(min(map_1), min(map_2)) - 1);
    axis_max = ceil(max(max(map_1), max(map_2)) + 1);
    
    xlim([axis_min axis_max])
    ylim([axis_min axis_max])
    hold on; refline(1)
    correlation = corrcoef(map_1, map_2);
    disp(strcat('Correlation coefficient: ', num2str(correlation(2))))
    
    figure; histogram(map_1, 'BinWidth', 0.5)
    hold on; histogram(map_2, 'BinWidth', 0.5)
    legend(map_1_label, map_2_label)
    