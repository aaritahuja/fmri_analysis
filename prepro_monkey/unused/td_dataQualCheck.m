function td_dataQualCheck(boldPath, runDirs, fprefix)

% Inputs:
% <boldPath> full path string of the upper level directory
% <runDirs> cell array of run folders
% <fprefix> string file prefix, usually f, af, raf, or wraf

% Based on rawQuant.m by CG
% Last modified 2/1/17 TMD


for ri = 1:length(runDirs)
    %fidx = ~cellfun(@isempty, strfind(subjFiles, sprintf('_run%d', runnums(ri))));
    
    rundir = fullfile(boldPath, runDirs{ri} );
    
    % Get the functional files for the current bold directory
    % May want to fix to pre-pend the files with f
    dirdata = dir(fullfile(rundir, [fprefix '*.nii']));
    
    % Convert struct to character array
    fnames = char(dirdata.name);
    % Make it the full file name
    dirnames = strcat(rundir, '/', fnames);
    
    % Three figures
    
    % tsdiffana
    % default position [0.3256    0.0084    0.3750    0.9105]
    hF(1) = figure('Units', 'normalized', 'Position', [0,0,.25,.9]);
    tsdiffana_sub(dirnames, 0, hF(1));
    
    % art_global
    hF(2) = figure('units', 'normalized', 'position', [.25, 0, .4, .8]);
    art_global_sub(dirnames, hF(2));
    
    % art_movie
    hF(3) = figure('units', 'normalized', 'position', [.65, 0, .35, .8]);
    art_movie_sub(1,dirnames,hF(3));
    
    fprintf('%s Run %d: %d files\n', boldPath, runDirs{ri}, length(fnames))
    inp = input('Next run & close figs (n), next run & open new figs (o), exit (x)? ', 's');
    switch inp
        case 'n'
            close(hF);
        case 'o'
            clear hF
        case 'x'
            break
        otherwise
            error('unrecognized option')
    end
end
   