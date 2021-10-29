function model_structure = model_all_conditions_subjects(MRI_directory, subjids, modeltypes, contrasts, ROI, varargin)

model_structure.struct_init = 1;
for s = 1:length(subjids)
    fprintf('Starting subject %d \n', subjids(s))
    %ROI_file = fullfile('/Users/Aaru/Documents/MRI_data', num2str(subjids(s)), 'models/localizer/higher_order_motion.mat');
    for m = 1:length(modeltypes)
        for c = 1:length(contrasts)
            subject = subjids(s);
            model = modeltypes{m};
            contrast = contrasts{c};
            
            if strcmp(ROI(end-7:end), 'resliced')
                ROI_file = fullfile(MRI_directory, num2str(subjids(s)), 'models/task_with_residuals', ROI);
            else
                ROI_file = fullfile(MRI_directory, 'RSA/ROIs', ROI);
            end
            
            if length(varargin) > 0
                subtraction_file = varargin{1};
                [beta, t, SE, CI, all_coords, translated_coords] = get_model_values(MRI_directory, subject, model, contrast, ROI_file, subtraction_file);
            else
                [beta, t, SE, CI, all_coords, translated_coords] = get_model_values(MRI_directory, subject, model, contrast, ROI_file);
            end
            struct_contrast = contrast(1:end-1);
            model_structure.(strcat('s', num2str(subject))).(model).(struct_contrast).beta = beta;
            model_structure.(strcat('s', num2str(subject))).(model).(struct_contrast).t = t;
            model_structure.(strcat('s', num2str(subject))).(model).(struct_contrast).SE = SE;
            model_structure.(strcat('s', num2str(subject))).(model).(struct_contrast).CI = CI;
            model_structure.(strcat('s', num2str(subject))).(model).(struct_contrast).coords = all_coords;
            model_structure.(strcat('s', num2str(subject))).(model).(struct_contrast).trans_coords = translated_coords;
            fprintf('%s contrast done \n', contrast)
        end
    end
    fprintf('Subject %d done \n', subject)
end

model_structure = rmfield(model_structure, 'struct_init');