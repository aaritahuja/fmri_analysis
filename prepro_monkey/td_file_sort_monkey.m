function td_file_sort_monkey(sessnum, subjectPath, anums, fnums, rnums)
% *** Assumes that the series number is in the folder name
% Inputs:
% <sessnum> session number, e.g. 2 - will be converted to '002'
% <subjectPath> full path (string) of the directory where the subject data
%   is stored, e.g. '/mnt/sdd1/mSEQ/subjects'
% <anums> array of series numbers of all the anatomicals
% <fnums> array of series numbers of all the functionals
% <rnums> array of series numbers of all the runs - will be numbered sequentially
% td_file_sort( ..., 'test') does not move or rename any files
% td_file_sort( ..., 'ftypes', <ftypes>) option to select which file types
%   should attempt to sort. Useful when restarting on a partially completed
%   session
% Example (where the first functional run, sequence 3, was aborted):
%   td_file_sort(101, '/mnt/sdd1/mSEQ/subjects', [1 2 10], 3:9, 4:9)

% Naming at the MRF seems to sometimes be just the subject number or the
%   study and then the subject number so have to compensate for that.c
% Change the wildcards to something like:
% files = dir('/datafiles2/SEQ1/101/bold/001/f*101-0015*');
% files = dir('/datafiles2/SEQ1/102/bold/001/f*102-0015*');

% Note, if the start of the file name is incorrectly named at the scanner,
%   can use something like td_fileRename to change everything before running
%   this


tic

basedir = subjectPath;
sessstr = sprintf('%.3d', sessnum);
subjPath = fullfile(basedir, sessstr);

% Rename subject folder from studyname_subjnum to just subjnum
finfo = dir(fullfile(basedir, ['*' int2str(sessnum) '*']));
if size(finfo,1) ~= 1
    error('Not 1 session dir found')
end
subjstr = finfo.name;

%%
% Before doing anything check one of the file names to make sure it
%   doesn't need renaming
tempsubjFiles = dir(fullfile(basedir, subjstr,'1*'));
%If rescanning and files start with a different number
tempsubjFiles = {tempsubjFiles.name};
% Used to take a particular folder, now just take the first one that fits
%   the naming criteria 
%idx = ~cellfun(@isempty, strfind(tempsubjFiles, sprintf('1_%.3d_', anums(1))));
tempaFiles = dir(fullfile(basedir, subjstr, tempsubjFiles{1}));
tempaFiles = {tempaFiles.name}; %contains . and ..
tempaFiles = tempaFiles(~ismember(tempaFiles, {'.' '..' '.DS_Store'}));

%%
% Do the actual renaming of the subject folder
[s,mess,messid] = movefile(fullfile(basedir, subjstr), subjPath);
if s ~= 1
    error(mess)
end
    
% Get names of directories
%   leave out . and .. by all of them starting with "1"
subjFiles = dir(fullfile(subjPath,'1*'));
%If not starting with 1
subjFiles = {subjFiles.name};

%% Move anatomicals
aPath = fullfile(subjPath, 'anatomy');
mkdir(aPath)

for i = 1:length(anums)
    idx = ~cellfun(@isempty, strfind(subjFiles, sprintf('1_%.3d_', anums(i))));

    [s,mess,messid] = movefile(fullfile(subjPath,subjFiles{idx}), aPath);
    if s ~= 1
        error(mess)
    end
end

%% Move & rename functionals
fPath = fullfile(subjPath, 'bold');
mkdir(fPath)
bad_path = fullfile(fPath, 'bad_runs');
mkdir(bad_path);

for i = 1:length(fnums)
    % Find the name of the sequence directory
    idx = ~cellfun(@isempty, strfind(subjFiles, sprintf('1_%.3d_', fnums(i))));
    
    if ~any(idx)
        error('no directory found')
    end
    
    % Just move it to the bold folder as-is
    %   Accounts for other files that may be in folder (eg. .img or .txt)
    [s,mess,messid] = movefile(fullfile(subjPath,subjFiles{idx}), fPath);
    if s ~= 1
        error(mess)
    end
    
    % If it's a run to be used, rename the folder name and the files within
    if ismember(fnums(i), rnums)
        % Figure out what run it is by the position in the array
        run = find(ismember(rnums, fnums(i)));
        
        seqPath = fullfile(fPath,subjFiles{idx});
        seqFiles = dir(fullfile(seqPath,'*.nii'));
        seqFiles = {seqFiles.name};
        
        % Rename the sequence folder
        rPath = fullfile(fPath, sprintf('%.3d', run));
        [s,mess,messid] = movefile(seqPath, rPath);
        if s ~= 1
            [s2,mess2,messid2] = movefile(seqPath, rPath);
            if s2 ~= 1
                error(mess)
            end
        end
        
        for k = 1:length(seqFiles)
            startstr = regexp(seqFiles{k}, subjstr, 'match');
            
            % Find the acquisition number (also takes the file ending)
            endstr = regexp(seqFiles{k},'_[0-9]{4}.nii$','match');
            
            % New file name (start with an f!)
            newfname = sprintf('f%ss%.3d_r%.3d%s', startstr{1}, fnums(i), ...
                run, endstr{1});
            
            
            [s,mess,messid] = movefile(fullfile(rPath, seqFiles{k}), ...
                fullfile(rPath, newfname) );
            if s ~= 1
                error(mess)
            end
        end
        fprintf('Moved %d files to %s\n', length(seqFiles), rPath)
    else
        movefile(fullfile(subjPath, 'bold', subjFiles{idx}), bad_path);
    end
end

%% Behavior
% Make the behavior folder and print a reminder to move files into it

fPath = fullfile(subjPath, 'behavior');
mkdir(fPath)
fprintf('Created %s REMINDER to copy behavioral files here\n', fPath)

%% Do the same for the data qual check figures folder
dPath = fullfile(subjPath, 'qualcheck');
mkdir(dPath)
fprintf('Created %s REMINDER to do and save data quality check figures here\n', dPath)

toc