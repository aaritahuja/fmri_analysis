function td_prepro(MRI_path, subjnum, varargin)

% function td_prepro(subjnum)
% Written for use with SPM12
% Input:
% <subjPath> full path to the "subjects" directory that should each contain
%   a numbered subject folder
% <subjnum> subject number. Assumes the folder name where the data are
%   stored are named with that number

% Example: 
% td_prepro(102)

% Output:
% The following files will be created in the upper level subject directory:
%   preprolog.txt - log file of times processes started and ended
%   batchslicetime.mat, batchrealign.mat, batchnorm.mat, batchsmooth.mat,
%   and batcht1norm.mat - all files that can be loaded into the SPM batch
%   processor
% SPM should also save spm_Date.ps that contains the graphs for realignment
% Update 12/11/17 to run on Oscar more generally

% NOTE: Some of the specific preprocessing options may have to be adjusted
% by study!

% Last modified 12/11/17 by Theresa M. Desrochers

% **********************************************************************
% Modify items in this section for individual subjects/studies

% Preprocessing options
% ---------------------
% IMPORTANT: Set which steps you would like to run for this subject.
% You can modify this with varargin 
prepro.slicetime.run = false; % Slice timing

prepro.realign.run = true; % Realignment
prepro.realign.runbyrun = true; % usually false

prepro.norm.run = true; % EPI normalization
prepro.norm.runbyrun = true; % usually false

prepro.smooth.run = true; % Smoothing
prepro.t1norm.run = true; % T1 normalization

% Default directories
% -------------------
% These will need to be modified/updated in general once per study
%subjnum = 101;
subjPath = fullfile(MRI_path, int2str(subjnum));
boldPath = fullfile(subjPath, 'bold');
rundirs = {'001' '002' '003' '004' '005' '006' '007' '008' '009' '010' '011' '012' '013'};

% Make sure spm is running
%startup_spm12
%spm fmri

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

% Set working directory to subject number so that all the graphs that are
%   automatically saved by SPM will go somewhere that makes sense
cd(subjPath)

% Open a log file (looks in the current directory)
if ~exist('preprolog.txt', 'file')
    fid = fopen('preprolog.txt','w');
    fprintf(fid, 'subjPath = %s\n', subjPath);
    fclose(fid);
end

% What's this?
% Warning: Run spm_jobman('initcfg'); beforehand 

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
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {'C:/Users/lab/Documents/MATLAB/spm12/tpm/TPM.nii'};
    %matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {'/Users/Aaru/Documents/MATLAB/spm12/tpm/TPM.nii'};
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
        error('not one t1foldname')
    end
    t1Path = fullfile(anatPath, t1foldname.name);
    t1fname = dir(fullfile(t1Path, '*t1_mprage*.nii'));
    if length(t1fname) ~= 1
        error('not one t1fname')
    end
    t1fPath = fullfile(t1Path, t1fname.name);
    
    % Subject's t1mprage is the one to normalize...
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = ...
        { t1fPath };

    % And the one to write
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = ...
        { t1fPath };
    
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
    %matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {'C:/Users/lab/Documents/MATLAB/spm12/tpm/TPM.nii'};
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {'/Users/Aaru/Documents/MATLAB/spm12/tpm/TPM.nii'};
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
        78 76 85];
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [2 2 2];
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;
    
    % Save the batch
    save(fullfile(subjPath, 'batcht1norm.mat'), 'matlabbatch')
    fprintf(fid, 'Saved matlabbatch to %s\n', ...
        fullfile(subjPath, 'batcht1norm.mat'));
    
    % Run the batch
    spm_jobman('run', matlabbatch);
    
    % Output complete
    fprintf(fid, 'T1 normalization completed on %s in %.3g min\n', ...
        t1fPath, toc/60);
    fprintf('T1 normalization completed on %s in %.3g min\n', ...
        t1fPath, toc/60)
    fclose(fid);
    clear matlabbatch
    
end %of t1 normaalization

%% Finishing

% Set the path up to the main data directory (outside of subjects)
cd('..')
cd('..')
fprintf('**********************\nDone preprocessing %d\n**********************\n', ...
    subjnum)