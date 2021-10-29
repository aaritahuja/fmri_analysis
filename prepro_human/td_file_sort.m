function td_file_sort(MRI_path, subjnum, anums, fnums, rnums)

%% Sort Converted .nii's after used MRIConvert on them
% *** Assumes that the series number is in the folder name
% Inputs:
% <subjnum> subject number, e.g. 101
% <subjectdir> full path (string) of the directory where the subject data
%   is stored, e.g. '/mnt/sdd1/mSEQ/subjects'
% <anums> array of series numbers of all the anatomicals
% <fnums> array of series numbers of all the functionals
% <rnums> array of series numbers of all the runs - will be numbered sequentially


%% Naming at the MRF seems to sometimes be just the subject number or the subeject number with preceding zeroes

% ******************Create five directories:
% anatomy
%   Move the anatomicals and don't change directory names for now
% bold
%   need an imput so can append run number to end of the directory
%   rename the epi files to start with f_experiment&subject_seqnum_some
%   descriptor of file name_scannumber
% behavior (empty for now)
%models (empty for now)
%preprocessing figures(empty for now)

tic

basedir = MRI_path;
subjPath = fullfile(basedir, int2str(subjnum));
finfo = dir(fullfile(basedir, ['*' int2str(subjnum) '*']));
subjstr = finfo.name;

% % Get names of directories
% %   leave out . and .. by all of them starting with "1"
subjFiles = dir(fullfile(subjPath,'1*'));
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

mprage_folder = dir(fullfile(aPath, '*mprage*'));
movefile(fullfile(aPath, mprage_folder.name), fullfile(aPath, 't1_mprage'));

%% Move & rename functionals
bPath = fullfile(subjPath, 'bold');
mkdir(bPath)
for i = 1:length(fnums)
    % Find the name of the sequence directory
    idx = ~cellfun(@isempty, strfind(subjFiles, sprintf('1_%.3d_', fnums(i))));
    
    % Just move it to the bold folder as-is
    %   Accounts for other files that may be in folder (eg. .img or .txt)
    [s,mess,messid] = movefile(fullfile(subjPath,subjFiles{idx}), bPath);
    if s ~= 1
        error(mess)
    end
    
    % If it's a run to be used, rename the folder name and the files within
    if ismember(fnums(i), rnums)
        % Figure out what run it is by the position in the array
        run = find(ismember(rnums, fnums(i)));
        
        seqPath = fullfile(bPath,subjFiles{idx});
        seqFiles = dir(fullfile(seqPath,'*.nii'));
        seqFiles = {seqFiles.name};
        
        % Rename the sequence folder
        rPath = fullfile(bPath, sprintf('%.3d', run));
        [s,mess,messid] = movefile(seqPath, rPath);
        if s ~= 1
            error(mess)
        end
        
        for k = 1:length(seqFiles)
            % Find the date
            subjprefix = strcat('0', subjstr);
            startstr = regexp(seqFiles{k}, subjprefix, 'match');
            
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
    end
end

%% Behavior
% Make the behavior folder and copy dgz files
bvPath = fullfile(subjPath, 'behavior');
mkdir(bvPath)
mkdir(fullfile(bvPath, 'localizer'))
mkdir(fullfile(bvPath, 'task'))
%Localizer files
bvlFiles = dir(fullfile('L:\projects\encounter\data\human', strcat('h', num2str(subjnum), '_motionloc*')));
bvlFiles = {bvlFiles.name};
for i = 1:length(bvlFiles)
    copyfile(fullfile('L:\projects\encounter\data\human', bvlFiles{i}), fullfile(bvPath, 'localizer'))
end
%Task files
task_files = {};
variants = {'native' 'counting' 'simulation' 'perception'};
for i = 1:length(variants)
    bvtFiles = dir(fullfile('L:\projects\encounter\data\human', strcat('h', num2str(subjnum), '_', variants{i}, '*')));
    bvtFiles = {bvtFiles.name};
    task_files = cat(2, task_files, bvtFiles);
end

for i = 1:length(task_files)
    copyfile(fullfile('L:\projects\encounter\data\human', task_files{i}), fullfile(bvPath, 'task'));
end

%% Models
% Make the models folder and subfolders
mPath = fullfile(subjPath, 'models');
mkdir(mPath)
mkdir (fullfile(mPath, 'localizer'))
mkdir (fullfile(mPath, 'task'))

%% Preprocessing figures
% Make the models folder and subfolders
pPath = fullfile(subjPath, 'preprocessing figures');
mkdir(pPath)
mkdir(pPath, 'before preprocessing')
mkdir (pPath, 'after realignment')
mkdir(pPath, 'after smoothing')

toc