function [sim_beta, per_beta, count_beta] = single_subject_similarity(MRI_directory, subjid, modeltype, ROI_file)

    [sim_beta, sim_SE, sim_CI] = get_model_values(MRI_directory, subjid, modeltype, 'simulation', ROI_file);
    [per_beta, per_SE, per_CI] = get_model_values(MRI_directory, subjid, modeltype, 'perception', ROI_file);
    [count_beta, count_SE, count_CI] = get_model_values(MRI_directory, subjid, modeltype, 'counting', ROI_file);
    
    sim_beta = rmmissing(sim_beta);
    per_beta = rmmissing(per_beta);
    count_beta = rmmissing(count_beta);

end

    
    
