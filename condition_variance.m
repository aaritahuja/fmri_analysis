function all_variances = condition_variance(MRI_directory, subjid, runs)

subject_directory = fullfile(MRI_directory, num2str(subjid));
bold_directory = fullfile(subject_directory, 'bold');
folder_number = '%03d';

all_variances = zeros(size(runs));

for i = 1:length(runs)
    %General opening settings
    run_directory = fullfile(bold_directory, sprintf(folder_number, runs(i)));
    files = dir(fullfile(run_directory, 'swraf*'));
    cd(run_directory)
    for imgi=1:length(files)
        vol = readnifti(files(imgi).name);
        % initialize
        if imgi==1
            fmridat = zeros([size(vol) length(files)]);
        end          
        fmridat(:,:,:,imgi) = vol;
    end
  run_variance = (var(fmridat, 0, 'all'))/100000;
  all_variances(i) = run_variance;
end
