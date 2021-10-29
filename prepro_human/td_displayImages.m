function td_displayImages(MRI_path, subjnum, runs, imgtype)

% Randomly picks two images of the type specificied <imgtype> from each of
%   the runs specified in <runs> for display. Program will wait for input 
%   to go to next image
% Image types:
% f - raw functionals
% a - slice timing corrected
% r - realigned
% w - normalized
% s - smoothed
% Example : td_displayImages(103, 1:6, 'f')

% 9/2/14 TMD

% Ensure that SPM 8 is on the path
%startup_spm12

spm fmri

% Adjust the file structure for your structure
%basedir = '/gpfs_home/dbasu1/data/pSEQ/subjects';

% basedir = 'Y:\data\pSEQ\subjects\'; % When using lab pc
%basedir = 'Z:\Human\pSEQ\subjects\'; % When using lab pc

subjPath = fullfile(MRI_path, int2str(subjnum));

for run = runs
    
    rPath = fullfile(subjPath, 'bold', sprintf('%.3d', run) );
    rFiles = dir(rPath);
    rFiles = {rFiles.name};
    
    tFiles = ~cellfun(@isempty, regexp(rFiles, sprintf('^%s*', imgtype) ));
    
    % From the files of the chosen type, choose two at random
    dispidxs = randsample(find(tFiles), 2);
    
    for i = 1:length(dispidxs)
        
        % Call the SPM code to display
        spm_image('Display', fullfile(rPath, rFiles{dispidxs(i)}))
        
        fprintf('Subj %d, run %d, img %d (%d of %d). ', subjnum, run, ...
            dispidxs(i), i, length(dispidxs))
        inp = input('Next image & overwrite figs (n) or exit (x)? ', 's');
        switch inp
            case 'n'
                % do nothing, just display the next
            case 'x'
                break
            otherwise
                error('unrecognized option')
        end
    end %of images to display
    
    if strcmpi(inp, 'x')
        break
    end
end %of runs
    