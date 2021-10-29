function model_structure = model_all_conditions_monkey(model_directory, contrasts, ROI)

model_structure.struct_init = 1;
%ROI_file = fullfile('/Users/Aaru/Documents/MRI_data', num2str(subjids(s)), 'models/localizer/higher_order_motion.mat');
%for m = 1:length(model_directories)
    for c = 1:length(contrasts)
        contrast = contrasts{c};
        %model_directory = model_directories{m};
        [beta, t, SE, CI, all_coords, translated_coords] = get_model_values_monkey(model_directory, contrast, ROI);
        
        struct_contrast = contrast(1:end-1);
        slashes = strfind(model_directory, '/');
        model_name = strcat('m', model_directory((slashes(end) + 1):end));
        
        model_structure.(model_name).(struct_contrast).beta = beta;
        model_structure.(model_name).(struct_contrast).t = t;
        model_structure.(model_name).(struct_contrast).SE = SE;
        model_structure.(model_name).(struct_contrast).CI = CI;
        model_structure.(model_name).(struct_contrast).coords = all_coords;
        model_structure.(model_name).(struct_contrast).trans_coords = translated_coords;
        fprintf('%s contrast done \n', contrast)
    end
%end

model_structure = rmfield(model_structure, 'struct_init');