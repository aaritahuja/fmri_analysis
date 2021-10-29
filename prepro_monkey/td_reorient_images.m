function td_reorient_images(fPath, fDirs, fprefix)

% <fPath> full string of the base path
% <fDirs> cell array of directories on fPath
% <fprefix> prefix of file. Typically 'f' for functional or 'walt' for
%   anatomical. Note that the flip direction is controled by this prefix

% Large swaths of code taken directly from spm_image.m case 'reorient'
% Copyright (C) 1994-2015 Wellcome Trust Centre for Neuroimaging
% John Ashburner

% Theresa M. Desrochers 2/22/17

% This is the userdata that would usually be entered for the
%   transformation in the display window. Will be used with spm_matrix. It
%   is a 12-element array. Note that "P" is used as a variable in
%   spm_image, so will use "B" here instead to stick with the naming
%   conventions in spm_image, rather than in spm_matrix.
% Return an affine transformation matrix
% FORMAT [A] = spm_matrix(P [,order])
% P(1)  - x translation (left-right)
% P(2)  - y translation (forward-back)
% P(3)  - z translation (up-down)
% P(4)  - x rotation about - {pitch} (radians)
% P(5)  - y rotation about - {roll}  (radians)
% P(6)  - z rotation about - {yaw}   (radians)
% P(7)  - x scaling
% P(8)  - y scaling
% P(9)  - z scaling
% P(10) - x affine
% P(11) - y affine
% P(12) - z affine

% This is set to rotate (pitch) by 90 degrees
%Original 
B = [-6 2 -6 1.7 0 0 1 1 1 0 0 0];
%B = [0 0 0 1.5708 0 0 1 1 1 0 0 0];
% Get the affine transformation matrix
M = spm_matrix(B);

if det(M) <= 0
    warning('This will flip the image(s)!')
end

% Save the matrix for reference
save(fullfile(fPath, [fprefix '_reorient.mat']), 'M');
fprintf('Saved affine transformation matrix to %s\n', ...
    fullfile(fPath, [fprefix '_reorient.mat']) )

for di = 1:length(fDirs)
    
    dirPath = fullfile(fPath, fDirs{di} );
    
    % Get the files for the current directory
    dirdata = dir(fullfile(dirPath, [fprefix '*.nii']));
    fnames = {dirdata.name};
    
    % For some reason fullfile was truncating the file name when doing it
    %   all at once in a cell array. cellfun didn't work either
    % There must be a limit to the number of characters that can be held in
    % a cell array, so do the adding of the path below! Looping didn't work
    % either

    spm_progress_bar('Init', length(fnames), ['Reorienting ' fDirs{di}],...
        'Images Complete');
    tic
    for fi = 1:length(fnames)
        
        % The reorienting machineary doesn't allow for renaming before the
        %   save, it does the orientation directly on the file. So save a
        %   copy of the file with a new name, and then do the operation on
        %   that copied file.
        
        % Prepend an 'l' to the file name
        flipPath = fullfile(dirPath, ['l' fnames{fi}]);
        % 'f' forces overwrite
        [status, message] = copyfile(fullfile(dirPath, fnames{fi}), flipPath, 'f');
        if status ~= 1
            error('copyfile of %s failed %s', flipPath, message)
        end
            
        % Initialize array (3rd dim is num of files);
        Mats = zeros(4,4,1);
        
        % Get the current voxel-to-world mappping of the image
        Mats(:,:,1) = spm_get_space(flipPath);
        
        % Set the new voxel-to-world mapping of the image (using M)
        spm_get_space(flipPath, M * Mats(:,:,1));
        
        
        spm_progress_bar('Set',fi);
    end
    spm_progress_bar('Clear');
    fprintf('Reoriented (''l'' prepend) %d %s images in %.2g min\n', ...
        length(fnames), fullfile(dirPath, [fprefix '*.nii']), toc/60 )
end
        