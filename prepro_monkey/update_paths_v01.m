function update_paths_v01(subjectnums, path_ccv, monkey, prepro)
%% Updating the paths from where the first level analysis occurs on /tmp using ccv. 
%% Changes path to the /gpfs mount and the file paths in our lab working directories on 
%% ccv so that data can then be processed, analyzed, etc 
% SPM has a built in function to do this: spm_changepath 

% v01: 09/18/2019 THM
% <subjectnums> Array of subject numbers- this comes from the array list in
% the ccv batch script and var name is: $SLURM_ARRAY_TASK_ID
% <path_ccv> variable in the shell script (.sh) is named $PATH_TMP - this
% is created on /tmp storage on ccv where we transfer our data and run our
% batch there. It also tacks on the $SLURM_JOB_ID to the filepath, so the
% filepath to files changes to: /tmp/$SLURM_JOB_ID = $PATH_TMP.


% Intialize SPM
spm defaults fmri
spm_jobman initcfg

%% Change these variables
%--------------------------

% Explanation
wraptxt = ('using spm changepath function to udpate file paths for /gpfs lab working dirs after running analysis on ccv /tmp space');

% This should be standard for a study, but want to double check
localpath = fullfile('/gpfs/data/tdesroch/monkey/', monkey);

%% Start the diary of command line output
% Because this is part of the batch script and save .out files, this info
% should be included there bc the function spits out text to the matlab
% command window. However, making a backup file just in case this info is
% needed for reference or paths get changed incorrectly, and then there are
% multiple records that can be double checked.

startDateTimeStr = datestr(now, 'yyyy-mm-dd-HH-MM');
diaryfname = sprintf('Diary_rms_updatepaths_%s', startDateTimeStr);

% Change grouppath here based on study name
diary(fullfile(localpath, diaryfname));

% Display current variables
disp(wraptxt)
fprintf('rms_updatepaths: %s\n', startDateTimeStr);
fprintf(['subjects %s\n' ...
    'groupAnalysisDir = %s\n'], ...
    mat2str(subjectnums), ...
    localpath);


tic

%% Constants
% ------------

for snum = subjectnums %subjectnums
    % Current location of the SPM.mat once it gets transferred back from
    % ccv stuff to our local group mount
    spmPath = fullfile('/gpfs/data/tdesroch/monkey/', monkey,...
         sprintf('%03d', 36), 'bold', strcat(prepro, '.mat'));
    % What we want the path to be
    newspmPath = fullfile('/gpfs/data/tdesroch/monkey', monkey);
    %fullfile('/gpfs/data/tdesroch/rmSEQ/subjects/', int2str(snum), '/'); % %This is what it was as of 12/16/2019, but it assumes the rest of the
    %path is /subjects/subjnum, so don't need that info here
    % Change the actual path
    spm_changepath(spmPath, path_ccv, newspmPath)
    %spm_changepath(spmPath, '/tmp/7764521/', newspmPath) % if we wrote out
    %path_ccv, it would look like this- tmp with a number associated with
    %the slurm job id
    
    fprintf('\nDone subject %s\n', num2str(snum))
end



fprintf('Done rms_updatepaths in %g min!\n', roundn(toc/60,-1))

diary off
end