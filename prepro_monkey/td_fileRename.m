function td_fileRename(sessPath, currPrefix, newPrefix)

% Written to specifically fix when the prefix to the file names is wrong
% <subjectPath> Full string path of the folder that contains all the scan
%   folders with scan files to be renamed
% <currPrefix> String of search for current prefix
% <newPrefix> String to replace currPrefix with.

% 5/5/17 TMD

% Find all the folders in the subject directory
dirinfo = dir(sessPath);
folders = {dirinfo.name};

for fi = 1:length(folders)
    
    if ismember(folders{fi}, {'.' '..' '.DS_Store' 'DICOM'})
        continue
    end
    
    folderPath = fullfile(sessPath, folders{fi});
    finfo = dir(folderPath);
    subjFiles = {finfo.name};
    
    for si = 1:length(subjFiles)
        
        if ismember(subjFiles{si}, {'.' '..'})
            continue
        end
        
        if regexp(subjFiles{si}, ['^' currPrefix])
            newfname = [newPrefix subjFiles{si}(length(currPrefix)+1:end)];
            
            srcPath = fullfile(folderPath, subjFiles{si});
            destPath = fullfile(folderPath, newfname);
            
            if strcmpi(srcPath, destPath)
                % Not all platforms are case sensitive, so work-around so
                %   to avoid "Cannot copy or move a file or directory onto
                %   itself" error.
                
                % Generate a temporary file name
                tname = tempname(folderPath);
                
                % Move to temporary name and then to final name
                movefile(srcPath, tname);
                [s,mess,messid] = movefile(tname, destPath);
            else
                [s,mess,messid] = movefile(srcPath, destPath);
            end
            
            if s ~= 1
                error(mess)
            end
        end
    end
    
    % Keep in mind that this count includes . and ..
    fprintf('Renamed %d files that started with %s to start with %s in %s\n', ...
        length(subjFiles), currPrefix, newPrefix, folderPath);
    
end