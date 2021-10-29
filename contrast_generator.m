function contrast_generator(SPM_path, positive_w1, negative_w, varargin)
%SPM_path: location of SPM file
%positive_w1: Positive vector
%negative_w: Negative vector
%varargin{1}: Optional second positive vector
%varargin{2}: Weight for optional second positive vector (if unspecified, a
%weight of 1 will be assigned)
    
%addpath('/shared/lab/projects/analysis/aarit/MATLAB')
    % Load SPM.mat
    load(fullfile(SPM_path, 'SPM.mat'));
    
    %Find all cells where experimental condition is 
    positive = double(contains(SPM.xX.name, positive_w1));
    
    if nargin >= 4
        positive_w2 = varargin{1};
        contrast_name = strcat(positive_w1, '+', positive_w2, '>', negative_w);
        if nargin == 5
            positive1 = double(contains(SPM.xX.name, positive_w2))*varargin{2};
        else
            positive1 = double(contains(SPM.xX.name, positive_w2));
        end
        negative_counter_weight = max(positive) + max(positive1);
        positive = positive + positive1;
    else
        negative_counter_weight = max(positive);
        contrast_name = strcat(positive_w1,'>', negative_w);
    end
          
    %Same as above but for control
    if strcmp(negative_w, '')
        negative = zeros(1, length(positive));
    else
        negative = double(contains(SPM.xX.name, negative_w));
        negative = negative*-1*negative_counter_weight;
    end
    
    if strcmp(negative_w, '')
        t_contrast = (positive/(sum(positive))) + (negative);
    else
        t_contrast = (positive/(sum(positive))) + (negative/(sum(negative))*-1);
    end
    disp(strcat("Contrast sum is ", num2str(sum(t_contrast))))
    %t_contrast = t_contrast/(sum(positive));
    
    %Roughly check it in the GUI, run manually
    matlabbatch{1}.spm.stats.con.spmmat = {fullfile(SPM_path, 'SPM.mat')};
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = contrast_name ;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = t_contrast;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.delete = 0;
    %spm_jobman('interactive',matlabbatch)
    spm_jobman('run',matlabbatch)
end