function zipscript(type, MRI_path, subjid, foldernames, sprefix, varargin)

%Function to zip intermediates, remove intermediates, or
%unzip intermediates.

% <type>
%If type=zip    all intermediate files (f, af, raf, wraf) are zipped
%   type=zipdel all intermediate files are deleted
%   type=unzip  all files are unzipped

% <foldernames> should be a cell array of string folder names
%   e.g. {'001' '002' '003' '004' '005'}
% <sprefix> is the string prefix that all file names will have and that the
%   original .nii files will start with. e.g. KAHN or SEQ1 or f (output of
%   dicomsort)
% zipscript(..., 'dir', <basedir>) gives the full path to the directory in which
%   to find the <foldernames>, usually /data/studyname/subjectnum/bold
%   otherwise will look for the folders in the current directory
% zipscript(..., 'prefix', <prefixes>) is a cell vector containing the
%   prefixes that should be operated on, e.g. {'ra'}

% Example call:

% zipscript('unzip', {'001' '002' '003' '004' '005'}, 'fSEQ', ...
%   'dir', '/gpfs_home/dbasu1/data/pSEQ/subjects/102/bold/', 'prefix', {'ra'});

% Example for zipping files
% zipscript('zip', {'001' '002' '003' '004' '005' '006'}, 'f') % Make
% correct bold directory is written below

%Written by Carolyn Ranti, 8.5.11
% Revised 11/14/13 Theresa Desrochers



% Set default prefixes
switch type
    case {'zip' 'unzip'}
        prefixes = { 'swra'};
         %prefixes = {'' 'a' 'ra' 'wra'};
        %prefixes = {'ra' 'wra'};
    case 'zipdel'
        %do *not* want to remove the original .nii files (f*)
        prefixes = {'a' 'ra' 'wra'};
    otherwise
        error('unrecognized type')
end

% basedir = cd;
%**************** CD to the bold data folder of current subject
%basedir= 'C:/Users/lab/Documents/MRI_data/6/bold';
%basedir = 'C:\Users\lab\Dropbox (Brown)\MRI_data\5\bold';
basedir= fullfile(MRI_path, subjid, 'bold') ;

argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'dir'
            argnum = argnum + 1;
            basedir = varargin{argnum};
        case 'prefix'
            argnum = argnum + 1;
            prefixes = varargin{argnum};
        otherwise
            error('unrecognized option')
    end
    argnum = argnum + 1;
end

for fi = 1:length(foldernames)
    fullfoldername = fullfile(basedir, foldernames{fi});
    cd(fullfoldername)
    if exist(fullfoldername, 'dir')
        switch type
            case 'zip'
                %zip files
                for i = 1:length(prefixes)
                    prefix = sprintf('%s%s*', prefixes{i}, sprefix);
                    %unixstring = sprintf('gzip %s', ...
                    %    fullfile(fullfoldername, prefix) );
                    %unix(unixstring);
                    zipfile = strcat('zip_', prefixes{i}, 'f');
                    zip(zipfile, prefix);
                    delete(prefix);
                    fprintf('done %s\n', fullfoldername, prefix);
                end
            case 'unzip'
                %unzip all files                
                for i = 1:length(prefixes)
                    prefix = sprintf('%s%s*', prefixes{i}, sprefix);
                    zipfile = strcat('zip_', prefixes{i}, 'f');
                    unzip(zipfile, fullfoldername);
                    fprintf('done %s\n', fullfoldername, prefix);
                end               
            otherwise
                error('zipscript:type', 'unrecognized type %s', type)
        end
    else
        warning('zipscript:folder', 'folder %s does not exist', ...
            fullfoldername);
    end
end

disp('Done!')
if strcmpi(type, 'zip')
    disp('Don''t forget to move spm_date.ps files.')
end