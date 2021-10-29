function td_prepro_rmsbatch(subjnum, path_ccv, path_subj, varargin)

% started from td_prepro(subjnum) in td_tools git folder. Copied to rmseq
% analysis folder to edit for batch processing on ccv
% Written for use with SPM12
% Input:
% <subjnum> subject number. Assumes the folder name where the data are
%   stored are named with that number
% <path_ccv> variable in the shell script (.sh) is named $PATH_TMP - this
% is created on /tmp storage on ccv where we transfer our data and run our
% batch there. It also tacks on the $SLURM_JOB_ID to the filepath, so the
% filepath to files changes to: /tmp/$SLURM_JOB_ID = $PATH_TMP.
% Can just use this variable to run the spm_changepath script to adjust the
% file paths back to our working directories on ccv for lab group
% <path_subj> variable in the shell script (.sh) is named $PATH_TMP_SUBJECTS - this
% is created on /tmp storage on ccv where we transfer our data and run our
% batch there. It also tacks on the $SLURM_JOB_ID to the filepath, so the
% filepath to files changes to: /tmp/$SLURM_JOB_ID/subjects = $PATH_TMP_SUBJECTS.

% Output:
% The following files will be created in the upper level subject directory:
%   preprolog.txt - log file of times processes started and ended
%   batchslicetime.mat, batchrealign.mat, batchnorm.mat, batchsmooth.mat,
%   and batcht1norm.mat - all files that can be loaded into the SPM batch
%   processor
% SPM should also save spm_Date.ps that contains the graphs for realignment

% Note that only editing the script right now to run slice timing
% correction because that is the step that is taking the longest to run on
% ccv

% Last modified 11/5/2019 by Theresa H. McKim

% **********************************************************************
% Modify items in this section for individual subjects/studies

% Preprocessing options
% ---------------------
% IMPORTANT: Set which steps you would like to run for this subject. 
prepro.slicetime.run = true; % Slice timing
prepro.realign.run = false; % Realignment
prepro.realign.runbyrun = false; % usually false
prepro.norm.run = false; %EPI normalization
prepro.norm.runbyrun = false; % usually false
prepro.smooth.run = false; % Smoothing
prepro.t1norm.run = false; % T1 normalization

% Default directories
% -------------------

%% Use line below if using path_ccv input and running with transfer to tmp
% space on ccv
paths.base = path_ccv;
% path to the subj dir
paths.subj = path_subj;

%% Avoid working directory errors
% Add the path where the scripts are so that can safely change directories
% to run the rest
% Use line below if running locally
%cd('/gpfs_home/tmckim/rmseq_analysis') % where github scripts are located

% Use next three lines below if using path_ccv input and running with transfer to tmp
% space on ccv
paths.scripts = fullfile(paths.base, 'scripts'); 
addpath(fullfile(paths.base, 'scripts')) %#update
cd(paths.scripts)

% For ccv:
subjPath = fullfile(paths.subj, int2str(subjnum));
boldPath = fullfile(subjPath, 'bold');
rundirs = {'001' '002' '003' '004' '005' '006'};


% These will need to be modified/updated in general once per study-
% previous
% subjPath = fullfile('/gpfs/data/tdesroch/rmSEQ/subjects', int2str(subjnum));
% boldPath = fullfile(subjPath, 'bold');
% rundirs = {'001' '002' '003' '004' '005' '006'};


% Make sure spm is running
% Intialize SPM
spm defaults fmri
spm_jobman initcfg
%spm_get_defaults('cmdline',true) % this ensures that if a dialog box pops up, it will show in the command line so you can see why something may fail in the .out file


argnum = 1;
while argnum  <= length(varargin)
    switch varargin{argnum}
        case 'opts'
            argnum = argnum + 1;
            prepro = varargin{argnum};
        otherwise
            error('unknown option')
    end
    argnum = argnum + 1;
end

% **********************************************************************

%% Removed this as of 11/6 bc want to run from where the scripts are located. Shouldn't be saving any graphs for the slice timing step
% Set working directory to subject number so that all the graphs that are
%   automatically saved by SPM will go somewhere that makes sense
% cd(subjPath)

% Open a log file (looks in the current directory)
if ~exist('preprolog.txt', 'file')
    fid = fopen('preprolog.txt','w');
    fprintf(fid, 'subjPath = %s\n', subjPath);
    fclose(fid);
end


%% Slice timing

if prepro.slicetime.run
    
    fid = fopen('preprolog.txt', 'a');
    fprintf(fid, 'Slice timing started at %s\n', datestr(now));
    fprintf('Slice timing started at %s\n', datestr(now))
    tic
    
    % Note, these parameters should be adjusted for each study! See
    %   FMRI_Prepro_SPM_...docx or the help under the Slice timing button for
    %   more details
    matlabbatch = {};
    %matlabbatch{1}.spm.temporal.st.scans = {'<UNDEFINED>'};
    matlabbatch{1}.spm.temporal.st.nslices = 38;
    matlabbatch{1}.spm.temporal.st.tr = 2;
    matlabbatch{1}.spm.temporal.st.ta = 1.94736842105263;
    matlabbatch{1}.spm.temporal.st.so = [2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37];
    matlabbatch{1}.spm.temporal.st.refslice = 1;
    matlabbatch{1}.spm.temporal.st.prefix = 'a';
    
    % Even though adding scans for all runs, have to be careful about adding
    %   recursively for all runs at once because there may be bold acquisition
    %   folders that want to leave out
    
    % Examples of how scans should be listed:
    %matlabbatch{1}.spm.temporal.st.scans = {'<UNDEFINED>'};
    % matlabbatch{1}.spm.temporal.st.scans = {
    %                                         {
    %                                         '/mnt/sdd1/mSEQ/subjects/101/bold/001/fmSEQ_101_20151106_s004_r001_0001.nii,1'
    %                                         '/mnt/sdd1/mSEQ/subjects/101/bold/001/fmSEQ_101_20151106_s004_r001_0002.nii,1'
    %                                         '/mnt/sdd1/mSEQ/subjects/101/bold/001/fmSEQ_101_20151106_s004_r001_0003.nii,1'
    %                                         '/mnt/sdd1/mSEQ/subjects/101/bold/001/fmSEQ_101_20151106_s004_r001_0004.nii,1'
    %                                         }
    %                                         }';
    scans = {};
    for ri = 1:length(rundirs)
        
        % Note that the file names must be the full path and must be a column
        %  vector
        runPath = fullfile(boldPath, rundirs{ri});
        
        % Find all the functional files for this run.
        % Assume that they are prefixed with an 'f'
        scanfnames = dir(fullfile(runPath, 'f*.nii'));
        
        % Need the full path, so concatinate the path and the name
        %   Note that fullfile does not add the last "/" so have to add that by
        %   hand
        scans = [scans; ...
            reshape(strcat(runPath, '/', {scanfnames.name}, ',1'), [], 1)];
    end
    
    % Add the scans back into the batch
    matlabbatch{1}.spm.temporal.st.scans = {scans}';
    
    % Save the batch
    save(fullfile(subjPath, 'batchslicetime.mat'), 'matlabbatch')
    fprintf(fid, 'Saved matlabbatch to %s\n', ...
        fullfile(subjPath, 'batchslicetime.mat'));
    
    % Run the batch
    spm_jobman('run', matlabbatch);
    
    % Output complete
    fprintf(fid, 'slice timing completed on %d files in %.3g min\n', ...
        length(scans), toc/60);
    fprintf('slice timing completed on %d files in %.3g min\n', ...
        length(scans), toc/60)
    fclose(fid);
    clear matlabbatch
end


%% Realignment

if prepro.realign.run
    
    fid = fopen('preprolog.txt', 'a');
    fprintf(fid, 'Realignment started at %s\n', datestr(now));
    fprintf('Realignment started at %s\n', datestr(now))
    tic
    
    % Note, these parameters should be adjusted for each study! See
    %   FMRI_Prepro_SPM_...docx or the help under the Realign,
    %   Realign: Estimate & Reslice button for more details
    matlabbatch = {};
%     matlabbatch{1}.spm.spatial.realign.estwrite.data = {
%         {
%         '/mnt/sdd1/mSEQ/subjects/101/bold/001/afmSEQ_101_20151106_s004_r001_0001.nii,1'
%         '/mnt/sdd1/mSEQ/subjects/101/bold/001/afmSEQ_101_20151106_s004_r001_0002.nii,1'
%         '/mnt/sdd1/mSEQ/subjects/101/bold/001/afmSEQ_101_20151106_s004_r001_0003.nii,1'
%         '/mnt/sdd1/mSEQ/subjects/101/bold/001/afmSEQ_101_20151106_s004_r001_0004.nii,1'
%         }
%         }';
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 3;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = {''};
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    
    % Even though adding scans for all runs, have to be careful about adding
    %   recursively for all runs at once because there may be bold acquisition
    %   folders that want to leave out
    
    data = {};
    nfiles = 0;
    for ri = 1:length(rundirs)
        
        % Note that the file names must be the full path and must be a column
        %  vector
        runPath = fullfile(boldPath, rundirs{ri});
        
        % Find all the functional files for this run.
        % Assume that they are prefixed with an 'af' from slice timing step
        scanfnames = dir(fullfile(runPath, 'af*.nii'));
        
        % Need the full path, so concatinate the path and the name
        %   Note that fullfile does not add the last "/" so have to add that by
        %   hand
        if prepro.realign.runbyrun
            data = reshape(strcat(runPath, '/', {scanfnames.name}, ',1'), [], 1);
            nfiles = nfiles + length(data);
            
            % Add the scans back into the batch
            matlabbatch{1}.spm.spatial.realign.estwrite.data = {data}';
            
            % Save the batch
            fname = fullfile(subjPath, sprintf('batchrealalign_%d.mat', ri));
            save(fname, 'matlabbatch')
            fprintf(fid, 'Saved matlabbatch to %s\n', fname);
            
            % Run the batch
            spm_jobman('run', matlabbatch);
            
            % Clear just the data, all the rest will stay the same
            matlabbatch{1}.spm.spatial.realign.estwrite.data = {};
        else
            data = [data; ...
                reshape(strcat(runPath, '/', {scanfnames.name}, ',1'), [], 1)];
        end
    end
    
    if ~prepro.realign.runbyrun
        % Add the scans back into the batch
        matlabbatch{1}.spm.spatial.realign.estwrite.data = {data}';
        
        % Save the batch
        save(fullfile(subjPath, 'batchrealalign.mat'), 'matlabbatch')
        fprintf(fid, 'Saved matlabbatch to %s\n', ...
            fullfile(subjPath, 'batchrealalign.mat'));
        
        nfiles = length(data);
        
        % Run the batch
        spm_jobman('run', matlabbatch);
    else
        fprintf(fid, 'RUN-BY-RUN\n');
        fprintf('RUN-BY-RUN\n');
    end
    
    
    %% Plotting realignment graph in batch mode- CCV
    %{
    %this is the dir where rp file is located- contains data for all 6 runs
    cd(fullfile(boldPath, '001'))
    
    rp = spm_load(spm_select('List','rp*.txt'));
    % SPM specific loading and selection of files
    % [files,dirs] = spm_select('List',direc,filt)
    %Return files matching the filter 'filt' and directories within 'direc'
    %  direc  - directory to search [Default: pwd]
    %  filt   - filter to select files with regexp, e.g. '^w.*\.img$' [Default: '.*'] 
    % files  - files matching 'filt' in directory 'direc'
    % dirs   - subdirectories of 'direc

     scaleme = [-3 3];

     printfig = figure;
     figure('Visible', 'off'); %CCV recommended to plot in batch mode
     set(printfig, 'Name', ['Motion parameters: subject ' num2str(subjnum) ]); % 'Visible', 'on');
     subplot(2,1,1);
     plot(rp(:,1:3));
     grid on;
     legend('x','y','z');
     ylim(scaleme);  % enable to always scale between fixed values as 
                     % set above
     title([ sprintf('Data from %s\n\n' , pwd) 'Translation (in mm)']);
           % 'interpreter', 'none');
     subplot(2,1,2);
     plot(rp(:,4:6)*180/pi);
     grid on;
     legend('x','y','z');
     ylim(scaleme);   % enable to always scale between fixed values as 
                      % set above
     title('Rotations (in dg)', 'interpreter', 'none');
     mydate = date;
     cd(subjPath)
     saveas(gcf, ['motionplot_' sprintf('%s', num2str(subjnum)) '_' mydate '.png']);
%}
     
    % Output complete
    fprintf(fid, 'Realignment completed on %d files in %.3g min\n', ...
        nfiles, toc/60);
    fprintf('Realignment completed on %d files in %.3g min\n', ...
        nfiles, toc/60)
    fclose(fid);
    clear matlabbatch
end

%% Normalization

if prepro.norm.run
    
    fid = fopen('preprolog.txt', 'a');
    fprintf(fid, 'Normalization started at %s\n', datestr(now));
    fprintf('Normalization started at %s\n', datestr(now))
    tic
    
    % Note, these parameters should be adjusted for each study! See
    %   FMRI_Prepro_SPM_...docx or the help under the Normalise,
    %   Normalise: Estimate & Reslice button for more details
    matlabbatch = {};
    if ~prepro.norm.runbyrun
        run1Path = fullfile(boldPath, rundirs{1});
        meanfname = dir(fullfile(run1Path, 'meanaf*.nii'));
        if length(meanfname) ~= 1
            error('more than one meanfname')
        end
        matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = ...
            { fullfile(run1Path, meanfname.name) };
    end
%     matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = {'/mnt/sdd1/mSEQ/subjects/101/bold/001/meanafmSEQ_101_20151106_s004_r001_0001.nii,1'};
%     matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = {
%         '/mnt/sdd1/mSEQ/subjects/101/bold/001/rafmSEQ_101_20151106_s004_r001_0001.nii,1'
%         '/mnt/sdd1/mSEQ/subjects/101/bold/001/rafmSEQ_101_20151106_s004_r001_0002.nii,1'
%         '/mnt/sdd1/mSEQ/subjects/101/bold/001/rafmSEQ_101_20151106_s004_r001_0003.nii,1'
%         };
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {'/gpfs/runtime/opt/spm/spm12/tpm/TPM.nii'};
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
        78 76 85];
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [2 2 2];
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;
    
    % Even though adding scans for all runs, have to be careful about adding
    %   recursively for all runs at once because there may be bold acquisition
    %   folders that want to leave out
    
    resample = {};
    nfiles = 0;
    for ri = 1:length(rundirs)
        
        % Note that the file names must be the full path and must be a column
        %  vector
        runPath = fullfile(boldPath, rundirs{ri});
        
        % Find all the functional files for this run.
        % Assume that they are prefixed with an 'raf' from realignment step
        scanfnames = dir(fullfile(runPath, 'raf*.nii'));
        
        % Need the full path, so concatinate the path and the name
        %   Note that fullfile does not add the last "/" so have to add that by
        %   hand
        if prepro.norm.runbyrun
            
            % Get the volume to compare to
            % Will only work if realigned run-by-run as well
            meanfname = dir(fullfile(runPath, 'meanaf*.nii'));
            if length(meanfname) ~= 1
                error('more than one meanfname')
            end
            matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = ...
                { fullfile(runPath, meanfname.name) };
            
            resample = reshape(strcat(runPath, '/', {scanfnames.name}, ',1'), [], 1);
            nfiles = nfiles + length(resample);
            
            % Add the scans back into the batch
            matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = resample;
            
            % Save the batch
            fname = fullfile(subjPath, sprintf('batchnorm_%d.mat', ri));
            save(fname, 'matlabbatch')
            fprintf(fid, 'Saved matlabbatch to %s\n', fname);
            
            % Run the batch
            spm_jobman('run', matlabbatch);
            
            % Clear just the vol & resample, all the rest will stay the same
            matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = {};
            matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = {};
        else
            resample = [resample; ...
                reshape(strcat(runPath, '/', {scanfnames.name}, ',1'), [], 1)];
        end
    end
    
    if ~prepro.norm.runbyrun
        % Add the scans back into the batch
        matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = resample;
        
        % Save the batch
        save(fullfile(subjPath, 'batchnorm.mat'), 'matlabbatch')
        fprintf(fid, 'Saved matlabbatch to %s\n', ...
            fullfile(subjPath, 'batchnorm.mat'));
        
        nfiles = length(resample);
        
        % Run the batch
        spm_jobman('run', matlabbatch);
    else
        fprintf(fid, 'RUN-BY-RUN\n');
        fprintf('RUN-BY-RUN\n');
    end
    
    % Output complete
    fprintf(fid, 'Normalization completed on %d files in %.3g min\n', ...
        nfiles, toc/60);
    fprintf('Normalization completed on %d files in %.3g min\n', ...
        nfiles, toc/60)
    fclose(fid);
    clear matlabbatch
end %of normalization

%% Smoothing
% Even if did realignment and normalization run-by-run, can do all runs at
%   once here because smoothing is independent of time

if prepro.smooth.run
    
    fid = fopen('preprolog.txt', 'a');
    fprintf(fid, 'Smoothing started at %s\n', datestr(now));
    fprintf('Smoothing started at %s\n', datestr(now))
    tic
    
    % Note, these parameters should be adjusted for each study! See
    %   FMRI_Prepro_SPM_...docx or the help under the Slice timing button for
    %   more details
    matlabbatch = {};
%     matlabbatch{1}.spm.spatial.smooth.data = {
%         '/mnt/sdd1/mSEQ/subjects/101/bold/001/wrafmSEQ_101_20151106_s004_r001_0001.nii,1'
%         '/mnt/sdd1/mSEQ/subjects/101/bold/001/wrafmSEQ_101_20151106_s004_r001_0002.nii,1'
%         '/mnt/sdd1/mSEQ/subjects/101/bold/001/wrafmSEQ_101_20151106_s004_r001_0003.nii,1'
%         };
    matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    
    % Even though adding scans for all runs, have to be careful about adding
    %   recursively for all runs at once because there may be bold acquisition
    %   folders that want to leave out
    
    data = {};
    for ri = 1:length(rundirs)
        
        % Note that the file names must be the full path and must be a column
        %  vector
        runPath = fullfile(boldPath, rundirs{ri});
        
        % Find all the functional files for this run.
        % Assume that they are prefixed with an 'f'
        scanfnames = dir(fullfile(runPath, 'wraf*.nii'));
        
        % Need the full path, so concatinate the path and the name
        %   Note that fullfile does not add the last "/" so have to add that by
        %   hand
        data = [data; ...
            reshape(strcat(runPath, '/', {scanfnames.name}, ',1'), [], 1)];
    end
    
    % Add the scans back into the batch
    matlabbatch{1}.spm.spatial.smooth.data = data;
    
    % Save the batch
    save(fullfile(subjPath, 'batchsmooth.mat'), 'matlabbatch')
    fprintf(fid, 'Saved matlabbatch to %s\n', ...
        fullfile(subjPath, 'batchsmooth.mat'));
    
    % Run the batch
    spm_jobman('run', matlabbatch);
    
    % Output complete
    fprintf(fid, 'Smoothing completed on %d files in %.3g min\n', ...
        length(data), toc/60);
    fprintf('Smoothing completed on %d files in %.3g min\n', ...
        length(data), toc/60)
    fclose(fid);
    clear matlabbatch
end

%% T1 Normalization

if prepro.t1norm.run
    
    fid = fopen('preprolog.txt', 'a');
    fprintf(fid, 'T1 normalization started at %s\n', datestr(now));
    fprintf('T1 normalization started at %s\n', datestr(now))
    tic
    
    % Note, these parameters should be adjusted for each study! See
    %   FMRI_Prepro_SPM_...docx or the help under the Normalise,
    %   Normalise: Estimate & Reslice button for more details
    matlabbatch = {};
    anatPath = fullfile(subjPath, 'anatomy');
    t1foldname = dir(fullfile(anatPath, '*t1_mprage*'));
    if length(t1foldname) ~= 1
        fprintf(fid,'more than one t1foldname');
        fprintf('more than one t1foldname')
    end
    
    % if they have more than one mprage and want to normalize both
    
    for i = 1:length(t1foldname)
        
        t1Path = fullfile(anatPath, t1foldname(i).name);
        t1fname = dir(fullfile(t1Path, '*t1_mprage*.nii'));
        %     if length(t1fname) ~= 1
        %         error('not one t1fname')
        %     end
        
        
        t1fPath = fullfile(t1Path, t1fname.name);
        
        % Subject's t1mprage is the one to normalize...
        matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = ...
            { t1fPath };
        
        % And the one to write
        matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = ...
            { t1fPath };
        
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {'/gpfs/runtime/opt/spm/spm12/tpm/TPM.nii'};
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
            78 76 85];
        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [2 2 2];
        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;
        
        % Save the batch
        fname = fullfile(subjPath, sprintf('batcht1norm_%d.mat', (i)));
        save(fname, 'matlabbatch')
        fprintf(fid, 'Saved matlabbatch to %s\n', fname);
            
        % Run the batch
        spm_jobman('run', matlabbatch);
        
        % Show how long it took for each mprage
        fprintf(fid, 'T1 normalization completed on %s in %.3g min\n', ...
            t1fPath, toc/60);
        fprintf('T1 normalization completed on %s in %.3g min\n', ...
            t1fPath, toc/60)
        
    end % for loop of multiple mprages
    
    
    % Output complete
    fprintf(fid, 'T1 normalization completed on %s in %.3g min\n', ...
        t1fPath, toc/60);
    fprintf('T1 normalization completed on %s in %.3g min\n', ...
        t1fPath, toc/60)
    fclose(fid);
    clear matlabbatch
    
end %of t1 normaalization

%% Finishing

% Set the path back to mSEQ
cd('/gpfs/data/tdesroch/rmSEQ')
fprintf('**********************\nDone preprocessing %d\n**********************\n', ...
    subjnum)