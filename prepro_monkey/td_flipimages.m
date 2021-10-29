function td_flipimages(fPath, fDirs, fprefix)

% <fPath> full string of the base path
% <fDirs> cell array of directories on fPath
% <fprefix> prefix of file. Typically 'f' for functional or 'walt' for
%   anatomical. Note that the flip direction is controled by this prefix
% Large swaths of code taken directly from imcalc.m
% a batch script to do imcalc transformations within SPM
% Copyright 2011, Robert J Ellis

% Theresa M. Desrochers 2/7/17

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
    %fpaths = fullfile(dirPath, {dirdata.name});
    
    % If wanted to use a mask, would need code such as this:
    % use_mask = 2;
    % maskch = 1;
    % Need to look up in the original code if want to impliment this
    %   functionality

    % Similarly don't need to change the data type
    % dtype_ch = 1;
    tic
    for fi = 1:length(fnames)
        
        % The format of this file name is a character array with the full path
        % of the file followed by ",1" - so need to create that
        % Note that this should be created here because it's too long to
        %   hold in a cell array
        file1 = fullfile(dirPath, [fnames{fi} ',1']);
        v1n = spm_vol(file1);
        v1  = spm_read_vols(spm_vol(file1));  % the actual volume
        v1(isnan(v1)) = 0;   % replace NaN if they are present
        
        
        % Original from imcalc.m:
        % vol = flipdim(v1,1);  % flip along the x-dimension
        
        % Updated for new matlab function and to flip in the y-dimension
        % (which is probably z given how the data were acquired)
        
        % Flip
        if strcmp(fprefix, 'f')
            % Functionals
            % dim 1 = flip L/R
            vol = flip(v1,3);
        else
            % Anatomical
            % flip along the first dimension, which is like z
            vol = flip(v1,1);
            % note the third dimension is a L/R flip, like x
        end
        
        % Prepend the new file with an "l" for "flip"
        v1n.('fname') = fullfile(dirPath, ['l' fnames{fi}]);
        
        v1n = spm_write_vol(v1n,vol);
        
    end
    fprintf('Flipped (''l'' prepend) %d %s images in %.2g min\n', ...
        length(fnames), fullfile(dirPath, [fprefix '*.nii']), toc/60 )
end
        