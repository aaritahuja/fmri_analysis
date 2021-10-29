function td_displayImages(boldPath, runs, imgtype)

% Randomly picks two images of the type specificied <imgtype> from each of
%   the runs specified in <runs> for display. Program will wait for input 
%   to go to next image
% <boldPath> string of the full 'bold' path
% <runs> row vector of runs to display images from. Note use numbers and
%   not folder names
% Image types:
% f - raw functionals
% a - slice timing corrected
% r - realigned
% w - normalized
% s - smoothed

% 2/2/17 TMD

% Ensure that SPM 12 is on the path
startup_spm12

for run = runs
    
    rPath = fullfile(boldPath, sprintf('%.3d', run) );
    rFiles = dir(rPath);
    rFiles = {rFiles.name};
    
    tFiles = ~cellfun(@isempty, regexp(rFiles, sprintf('^%s*', imgtype) ));
    
    % From the files of the chosen type, choose two at random
    dispidxs = randsample(find(tFiles), 2);
    
    for i = 1:length(dispidxs)
        
        % Call the SPM code to display
        spm_image('Display', fullfile(rPath, rFiles{dispidxs(i)}))
        
        fprintf('%s, run %d, img %d (%d of %d). ', boldPath, run, ...
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
    