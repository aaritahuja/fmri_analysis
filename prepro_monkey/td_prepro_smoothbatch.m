function td_prepro_smoothbatch(path_ccv, sessPath, anatPath, varargin)
%12/22/2019 NYR
% Changing based on td_prepro_rmsbatch by THM to be able to run as batches
% through CCV for faster reading and writing of output files.
% All of the batch processes will call this function, but there will be an 
% individual .sh file to run them 
% function td_prepro(subjnum)
% Written for use with SPM12
% Input:
% <prepro> Struct with the following fields to indicate whether an analysis
%   should run. They default to false.
% [required] prepro.t1folder to be the folder name *beyond* the
%   /anatomy/ folder where the t1 to be used is stored.
%     prepro.slicetime.run = false; % Slice timing
%     prepro.realign.run = false; % Realignment
%     prepro.realign.runbyrun = false; % usually false
%     prepro.t1norm.run = false; % T1 normalization
%       [required] prepro.t1norm.fname = just the filename, assumed to be
%       within prepro.t1folder
%     prepro.epicoreg.run = false; %coregister EPIs to native GM
%     prepro.epicoreg.runbyrun = false;
%     prepro.epinorm.run = false; % EPI normalization
%     prepro.epinorm.runbyrun = false; % usually false
%     prepro.smooth.run = false; % Smoothing
% <sessPath> Full string path containing data for a particular
%   subject/session. Assumes functional data are in <sessPath>/bold
% <runDirs> Numerical array of sub-directories in <sessPath>/bold to be
%   processed. Makes the assumption that run directories have all been
%   renamed to be 3 digits preceeded by filling zeros, e.g. 001.
% Output:
% The following files will be created in the upper level subject directory:
%   preprolog.txt - log file of times processes started and ended
%   batchslicetime.mat, batchrealign.mat, batchnorm.mat, batchsmooth.mat,
%   and batcht1norm.mat - all files that can be loaded into the SPM batch
%   processor
% SPM should also save spm_Date.ps that contains the graphs for realignment

% Last modified 05/22/17 by Theresa M. Desrochers

% **********************************************************************
%Make prepro field exist with mean mprage folder
prepro.t1folder = 't1_mprage_mean/';

% Default Preprocessing options
% -----------------------------
if ~isfield(prepro, 'slicetime')
    prepro.slicetime.run = false; % Slice timing
end
if ~isfield(prepro, 'realign')
    %Set to true so it defaults to realignment
    prepro.realign.run = false; % Realignment
    prepro.realign.runbyrun = false; % usually false
else
    if ~isfield(prepro.realign, 'runbyrun')
        prepro.realign.runbyrun = false;
    end
end
if ~isfield(prepro, 't1norm')
    prepro.t1norm.run = false; % T1 normalization
end
if ~isfield(prepro, 'epicoreg')
    prepro.epicoreg.run = false; %coregister EPIs to native GM
    prepro.epicoreg.runbyrun = false;
else
    if ~isfield(prepro.epicoreg, 'runbyrun')
        prepro.epicoreg.runbyrun = false;
    end
end
if ~isfield(prepro, 'epinorm')
    prepro.epinorm.run = false; % EPI normalization
    prepro.epinorm.runbyrun = false; % usually false
else
    if ~isfield(prepro.epinorm, 'runbyrun')
        prepro.epinorm.runbyrun = false;
    end
end
if ~isfield(prepro, 'smooth')
    prepro.smooth.run = true; % Smoothing
end

% Make sure spm is running
% startup_spm12
% spm fmri
%% Use line below if using path_ccv input and running with transfer to tmp
% space on ccv
paths.base = path_ccv;
% path to the subj dir
paths.sess = sessPath;

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
boldPath = fullfile(sessPath, 'bold');
%runDirNums = [1:12];
% We need to get the amounts of dirs directly from the folders
% grab dir names
folders = dir(boldPath);
% Use regexp to get number of folders that start with 0
sessDirs = sum(cell2mat(regexp({folders.name},'^0\w*')));
runDirNums = 1:sessDirs;


argnum = 1;
while argnum  <= length(varargin)
    switch varargin{argnum}
        %case 'opts'
            %argnum = argnum + 1;
            %prepro = varargin{argnum};
        otherwise
            error('unknown option')
    end
    argnum = argnum + 1;
end

% **********************************************************************

%% Set Paths
%  ---------

% Set working directory to subject number so that all the graphs that are
%   automatically saved by SPM will go somewhere that makes sense
%origWorkPath = cd;
%cd(sessPath)
%boldPath = fullfile(sessPath, 'bold');
if ismac
    atlasPath = '/Volumes/BM_DesrochersLab/Monkey/atlasfiles';
    templatePath = fullfile(atlasPath, '112RM-SL_T1.nii');
elseif isunix
    atlasPath = '/gpfs/data/tdesroch/monkey/atlasfiles';
    templatePath = fullfile(atlasPath, '112RM-SL_T1.nii');
else
    atlasPath = 'Z:\Monkey\atlasfiles';
    templatePath = fullfile(atlasPath, '112RM-SL_T1.nii');
end

%Had to add directly in batch script 
%anatPath = fullfile(sessPath, 'anatomy');

% Used to assume that there was only on directory that contained a
%   t1mprage, but that's not the case for averaging, so now expect it as an
%   input.
% t1foldname = dir(fullfile(anatPath, '*t1_mprage*'));
% if length(t1foldname) ~= 1
%     %error('not one t1foldname')
%     warning('not one t1foldname')
%     t1foldname = t1foldname(1);
% end

% Check to make sure that the t1 path makes sense
T1DirPath = fullfile(anatPath, prepro.t1folder);
if ~exist(T1DirPath, 'dir')
    error('Full T1DirPath is invalid. Check contents of input prepro.t1folder')
end

% Open a log file (looks in the current directory)
if ~exist('preprolog.txt', 'file')
    fid = fopen('preprolog.txt','w');
    fprintf(fid, 'sessPath = %s\n', sessPath);
    fclose(fid);
end

% What's this?
% Warning: Run spm_jobman('initcfg'); beforehand 

%% Slice timing
%{
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
    for ri = 1:length(runDirs)
        
        % Note that the file names must be the full path and must be a column
        %  vector
        runPath = fullfile(boldPath, runDirs{ri});
        
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
    save(fullfile(sessPath, 'batchslicetime.mat'), 'matlabbatch')
    fprintf(fid, 'Saved matlabbatch to %s\n', ...
        fullfile(sessPath, 'batchslicetime.mat'));
    
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
%}

%% Realignment

if prepro.realign.run
    
    fid = fopen('preprolog.txt', 'a');
    fprintf(fid, 'Realignment started at %s\n', datestr(now));
    fprintf('Realignment started at %s\n', datestr(now))
    tic
    
    % Set up batch
    % ------------
    % See the help under the Realign,
    %   Realign: Estimate & Reslice button for more details
    matlabbatch = {};
    
    % Estimation Options
    % ------------------
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 1;
    
    % Separation
    %   mm between the points sampled in the reference image
    %   Load an example image to get the dimensions - this is not done run
    %   by run!
    runDir1 = sprintf('%.3d', runDirNums(1));
    scanfnames = dir(fullfile(boldPath, runDir1, 'lf*001.nii'));
    v = spm_vol(fullfile(boldPath, runDir1, scanfnames.name));
    p = spm_imatrix(v.mat);
    % Parameters 7:9 are the x,y,z distances, see spm_matrix for help
    dims = abs(roundn(p(7:9), -1));
    if any(diff(dims))
        error('dimensions not the same')
    end
    % Take the first dimension as the separation
    % Must be 1-by-1 real number array
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = dims(1);
    
    % Smoothing
    %   Full Width Half Max of the Gaussian smoothing kernel (mm) applied
    %   to the images before estmating the realignment parameters
    %   Use 5 for humans, try 2 for monkeys
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 2;
    
    % Num Passes
    %   Register to mean: A two pass procedure is used in order to register
    %   the images to the mean of the images after the first realignment
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    
    % Interpolation: 2nd Degree B-spline
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    
    % Wrapping - currently set for none, may want to change
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    
    % Weighing - none
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = {''};
    
    % Reslice (write) Options
    % -----------------------
    
    % Resliced images: All Images + Mean Image
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    
    % Interpolation: 4th Degree B-spline
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    
    % Set up batch data
    
    % Even though adding scans for all runs, have to be careful about adding
    %   recursively for all runs at once because there may be bold acquisition
    %   folders that want to leave out
    
    %     matlabbatch{1}.spm.spatial.realign.estwrite.data = {
%         {
%         '/mnt/sdd1/mSEQ/subjects/101/bold/001/afmSEQ_101_20151106_s004_r001_0001.nii,1'
%         '/mnt/sdd1/mSEQ/subjects/101/bold/001/afmSEQ_101_20151106_s004_r001_0002.nii,1'
%         '/mnt/sdd1/mSEQ/subjects/101/bold/001/afmSEQ_101_20151106_s004_r001_0003.nii,1'
%         '/mnt/sdd1/mSEQ/subjects/101/bold/001/afmSEQ_101_20151106_s004_r001_0004.nii,1'
%         }
%         }';
    
    data = {};
    nfiles = 0;
    for ri = 1:length(runDirNums)
        
        % Note that the file names must be the full path and must be a column
        %  vector
        runPath = fullfile(boldPath, sprintf('%.3d', runDirNums(ri)));
        
        % Find all the functional files for this run.
        % Assume that they are prefixed with an 'af' from orienting step
        scanfnames = dir(fullfile(runPath, 'lf*.nii'));
        
        % Need the full path, so concatinate the path and the name
        %   Note that fullfile does not add the last "/" so have to add that by
        %   hand
        if prepro.realign.runbyrun
            data = reshape(strcat(runPath, filesep, {scanfnames.name}, ',1'), [], 1);
            nfiles = nfiles + length(data);
            
            % Add the scans back into the batch
            matlabbatch{1}.spm.spatial.realign.estwrite.data = {data}';
            
            % Save the batch
            fname = fullfile(sessPath, sprintf('batchrealign_%d.mat', ri));
            save(fname, 'matlabbatch')
            fprintf(fid, 'Saved matlabbatch to %s\n', fname);
            
            % Run the batch
            spm_jobman('run', matlabbatch);
            
            % Clear just the data, all the rest will stay the same
            matlabbatch{1}.spm.spatial.realign.estwrite.data = {};
        else
            data = [data; ...
                reshape(strcat(runPath, filesep, {scanfnames.name}, ',1'), [], 1)];
        end
    end
    
    if ~prepro.realign.runbyrun
        % Add the scans back into the batch
        matlabbatch{1}.spm.spatial.realign.estwrite.data = {data}';
        
        % Save the batch
        save(fullfile(boldPath, 'batchrealign.mat'), 'matlabbatch')
        fprintf(fid, 'Saved matlabbatch to %s\n', ...
            fullfile(boldPath, 'batchrealign.mat'));
        
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

%% T1 Normalization

if prepro.t1norm.run
    
    % This step comes next according to aamod_SegCoregNorm2d.m
    % Note that I am departing from using matlabbatch here in the interest
    %   of using the same options that Danny Mitchell uses
    
    fid = fopen('preprolog.txt', 'a');
    fprintf(fid, 'T1 normalization started at %s\n', datestr(now));
    fprintf('T1 normalization started at %s\n', datestr(now))
    tic
       
    % See the help under the Normalise,
    %   Normalise: Estimate & Reslice button for some more details
    %   Note that not using the batch options here, but calling the spm
    %   funcitons directly
    
    % Set Paths
    % ---------
    
    % Get the name of the t1 to use and don't assume anything about its name 
    t1fname = dir(fullfile(T1DirPath, prepro.t1norm.fname));
    if length(t1fname) ~= 1
        error('not one t1fname')
    end
    % t1fname.name ends in .nii
    T1name = t1fname.name;
    T1Path = fullfile(T1DirPath, T1name);
    mT1name = sprintf('m%s',T1name);
    mT1Path = fullfile(T1DirPath, mT1name);
    mT1Path_bet = strrep(mT1Path,'.nii','_betT1.nii');
    
    %% 1) Bias-corrected T1
    
    % Use the normalise functions to remove the intesity gradients
    % Makes the image look nicer, and may aid the initial skull stripping
    % Adds an 'm' in front of the image name
    
    % Estimate options
    % ----------------
    
    %defs.estimate.weight = '';
    estopts.tpm=char({...
        fullfile(atlasPath,'gm_priors_ohsu+uw.nii,1'),...
        fullfile(atlasPath,'wm_priors_ohsu+uw.nii,1'),...
        fullfile(atlasPath,'csf_priors_ohsu+uw.nii,1'),...
        });
    estopts.ngaus = [2 2 2 8]; % default = [2 2 2 4]
    estopts.biasfwhm= 30; % default = 60-75; reduce for monkey brain?
    estopts.samp = 2; % default=3; McLaren 2010 used 2
    estopts.warpreg = 1; % default =  1
    % warpco can stay large here for speed
    estopts.warpco = 25; % default = 25; reduce later for monkey brain?
    estopts.biasreg = 0.0001; % default = 0.0001
    estopts.msk='';
    
    % Write options
    % -------------
    
    % only write out attenuation corrected image
    writeopts.biascor = 1;
    writeopts.GM  = [0 0 0];
    writeopts.WM  = [0 0 0];
    writeopts.CSF = [0 0 0];
    writeopts.cleanup = 0;
    estopts.regtype=''; % turn off affine
    
    if ~exist(mT1Path, 'file')
        fprintf('\nCreating bias-corrected T1...')
        warning off all
        out = spm_preproc(T1Path, estopts);
        warning on all
        [sn] = spm_prep2sn(out); % don't write deformation files
        spm_preproc_write(sn, writeopts);
    else
        fprintf('\nBias-corrected T1 found.')
    end
    
    %% 2) Rough skull strip w/ bet
    
    % If you start too far away from the template when normalizing, then it
    %   fails. This rough skull strip is to help with co-registration in the
    %   next step. It uses FSL.
    if ~exist(mT1Path, 'file')

        fprintf('\nSkull-stripping structural using BET...')
    
    % f = fractional intensity threshold; smaller values give larger brain
    %   outline estimates
    % g = threshold gradient; positive values give larger brain outline at
    %   bottom, smaller at top
    % Original call didn't include any flags except for -v, so probably
    %   called with defaults: f=0.5 and g=0
        if ispc
            %If running on local machine use this path
            cmd = sprintf('/usr/local/fsl/bin/bet %s %s  -f 0.4 -g 0 -v', mT1Path, mT1Path_bet);
            call_fsl(cmd); %FSL matlab command
        elseif isunix
            %If running it on unix like would be the case for CCV, use this
            %instead
            cmd = sprintf('/gpfs/runtime/opt/fsl/6.0.0/bin/bet %s %s  -f 0.4 -g 0 -v', mT1Path, mT1Path_bet);
            call_fsl(cmd); %FSL matlab command
        end
    % Unzip everything
    %   r = recursive (if want to use directory, this case just 1 file)
    %   f = force overwrite of output file
        cmd = sprintf('gunzip -f %s.gz', mT1Path_bet);
        unix(cmd);
    else 
        fprintf('\nSkull stripped BET found.')
    end
    
    %% 3) Coregister mT1Path_bet to template
    
    % This will use the bias-corrected, skull stripped (bet) image to
    %   coregister to the template image (that also has no skull). This
    %   is necessary to get the images in roughly the same space (without
    %   warping) so that the normalization has a good place to start.
    
    % Defaults
    % --------
    % This does use SPM's batch functions
    coregjob{1}.spatial{1}.coreg{1}.estimate.other = {''};
    coregjob{1}.spatial{1}.coreg{1}.estimate.eoptions.cost_fun = 'nmi';
    coregjob{1}.spatial{1}.coreg{1}.estimate.eoptions.sep = [4 2 1]; % default is [4 2]
    coregjob{1}.spatial{1}.coreg{1}.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001]./2; % half of the defaults
    coregjob{1}.spatial{1}.coreg{1}.estimate.eoptions.fwhm = [7 7];
    % This option appears not to exist anymore
    %coregjob{1}.spatial{1}.coreg{1}.estimate.eoptions.graphics=0;
    
    % Coregister to template
    coregjob{1}.spatial{1}.coreg{1}.estimate.ref = {templatePath};
    
    % The thing to move around, the bias-corrected skull stripped image
    %   Could use just the bias-corrected (mT1Path), but Mitchell said it
    %   doesn't work as well
    coregjob{1}.spatial{1}.coreg{1}.estimate.source = {mT1Path_bet};
    
    % Other images to translate
    % Original script generates a list. Not needed here
    %   strucs=cellstr(spm_select('FPlist',fullfile(subjpath,'structurals'),...
    %   sprintf('^m?%s.*.nii$',aap.directory_conventions.rawdataafterconversionprefix)));
    % Coregister the original bias-corrected (not skull stripped)
    coregjob{1}.spatial{1}.coreg{1}.estimate.other = {mT1Path};
    
    % Run the SPM job
    % ---------------
    fprintf('\nCoregistering structural images to template...')
    spm_jobman('run', coregjob)
    
    %% 4) Normalise coregistered mT1Path to McLaren
    
    % Uses the bias-corrected, co-registered image so it's relatively
    %   speaking starting in the same space
    
    % Normalize & segment options (norm to McLaren in SegCoregNorm2d)
    
    % Estimate options
    % ----------------
    estopts.biasreg = 0.0001; % default = 0.0001
    estopts.biasfwhm= 30; % default = 60-75; reduce for monkey brain?
    estopts.tpm=char({...
        fullfile(atlasPath,'gm_priors_ohsu+uw.nii,1'),...
        fullfile(atlasPath,'wm_priors_ohsu+uw.nii,1'),...
        fullfile(atlasPath,'csf_priors_ohsu+uw.nii,1'),...
        });
    estopts.ngaus = [2 2 2 8]; % default = [2 2 2 4]
    estopts.msk='';
    estopts.regtype='subj';    % turn on affine: subj corresponds to "average sized template" as used by McLaren 2010
    estopts.samp = 2; % default=3; McLaren 2010 used 2
    estopts.warpco = 12; % default = 25; reduce for monkey brain?
    estopts.warpreg = 1; % default =  1
    
    % Do the Estimation
    % -----------------
    
    % Remember, this is on the bias-corrected, now coregistered image
    out = spm_preproc(mT1Path, estopts);
    [sn,isn] = spm_prep2sn(out); % don't write deformation files yet...
    clear out
    
    % Write options
    % -------------
    
    % Write native and modulated segmentations
    %   Bias correction already saved, so don't write it here
    writeopts.biascor = 0; 
    writeopts.GM  = [1 0 1];
    writeopts.WM  = [1 0 1];
    writeopts.CSF = [1 0 1];
    writeopts.cleanup = 1; % "light"
    
    % Do the Write
    % ------------
    spm_preproc_write(sn, writeopts);
    
    % Write forward and inverse deformation files
    % -----
    
    snmatname = [spm_str_manip(mT1Path,'sd') '_seg_sn.mat'];
    savefields(snmatname, sn);
    savefields([spm_str_manip(mT1Path,'sd') '_seg_inv_sn.mat'], isn);
    
    %% 5) Make maximum-liklihood tissue maps and brain mask
    % aamod_SegCoregNorm2d.m (moved from aamod_makeROIs)
    
    fprintf('\nClassifying tissues...')
    
    % T1name ends in .nii
    gm = spm_select('FPlist', T1DirPath, sprintf('^c1.*%s.*.nii$',T1name(1:end-4)));
    wm = spm_select('FPlist', T1DirPath, sprintf('^c2.*%s.*.nii$',T1name(1:end-4)));
    csf = spm_select('FPlist', T1DirPath, sprintf('^c3.*%s.*.nii$',T1name(1:end-4)));
    
    Vg=spm_vol(gm); 
    Yg=spm_read_vols(Vg); 
    if isfield(Vg,'pinfo'),
        Vg=rmfield(Vg,'pinfo'); 
    end
    Vw=spm_vol(wm); 
    Yw=spm_read_vols(Vw); 
    if isfield(Vw,'pinfo'),
        Vw=rmfield(Vw,'pinfo'); 
    end
    Vc=spm_vol(csf); 
    Yc=spm_read_vols(Vc);
    if isfield(Vc,'pinfo'),
        Vc=rmfield(Vc,'pinfo'); 
    end
    
    Vg.fname = fullfile(T1DirPath,'GM.nii');
    spm_write_vol(Vg, Yg>Yw & Yg>Yc & Yg>0.05);
    Vw.fname = fullfile(T1DirPath,'WM.nii');
    spm_write_vol(Vw, Yw>Yg & Yw>Yc & Yw>0.05);
    Vc.fname = fullfile(T1DirPath,'CSF.nii');
    spm_write_vol(Vc, Yc>Yw & Yc>Yg & Yc>0.05);
    clear Yg Yw Yc
        
    % Make brain mask
    % ---------------
    segmentations = spm_select('FPlist', T1DirPath, sprintf('^c\\d.*%s.*.nii$',T1name(1:end-4)));
    V = spm_vol(segmentations);
    Y = spm_read_vols(V);
    if isfield(V,'pinfo'),
        V = rmfield(V,'pinfo'); 
    end
    
    % Note that I'm saving to the anatomy path and not the T1path
    V(1).fname = fullfile(anatPath,'ssbrainmask.nii');
    if isfield(V,'pinfo'); 
        V = rmfield(V,'pinfo'); 
    end
    brainmask = sum(Y,4)>0.5;
    brainmask = spm_dilate(double(brainmask));
    spm_write_vol(V(1), brainmask);
        
    %% 6) Write normalised T1 and segmentations
    
    % Write options
    %   gathered from other places in aamod_SegCoregNorm2d.m
    % From aa3.userscript_PPcontrol_120314.m
    %aap.spm.defaults.normalise.write.vox=[2 2 2]; % resolution of normalised EPIs
    %aap.spm.defaults.normalise.write.bb=[-35 -30 -10; 35 55 45];
    %defs = aap.spm.defaults.normalise;
    %defs.write.bb = aap.spm.defaults.normalise.write.bb;
    defs.write.bb = [-35 -30 -10; 35 55 45];
    defs.write.vox = [0.5 0.5 0.5];
    
    spm_write_sn(char(segmentations, mT1Path, Vg.fname, Vw.fname, Vc.fname, V.fname), ...
        snmatname, defs.write);
    
    %% 7) Skull-strip T1 
    % (nicer viewing without hyperintensity at ears)
    V = spm_vol(mT1Path); 
    Y = spm_read_vols(V);
    if isfield(V,'pinfo'),
        V=rmfield(V,'pinfo'); 
    end
    
    % Note that I'm saving to the anatomy path and not the T1path
    V.fname = fullfile(anatPath, 'ssT1.nii');
    spm_write_vol(V, Y.*brainmask);
    clear Y brainmask
    
    % Write the normalized, skull stripped T1 also to the anatomy folder
    %   creates wssT1.nii
    defs.write.bb = [-35 -30 -10; 35 55 45];
    defs.write.vox = [0.5 0.5 0.5];
    spm_write_sn(fullfile(anatPath, 'ssT1.nii'), snmatname, defs.write);
    
    % Output complete
    fprintf(fid, 'T1 normalization completed on %s in %.3g min\n', ...
        T1Path, toc/60);
    fprintf('T1 normalization completed on %s in %.3g min\n', ...
        T1Path, toc/60)
    fclose(fid);
    %clear matlabbatch
    
end %of t1 normalization

%% EPI Coregistration

% In t1 normalization, recall that the "native" images were coregistered to
%   the template. In order to be able to apply the same normalization to
%   the EPIs as was applied to the anatomical, it is necessary to start
%   them off in the same place - coregistered to the template. However,
%   rather than coregistering to the template, it will be easier to
%   coregister to the native gray matter segmentation.
% Taken from "coregister epis to native GM segmentation" of
%   aamod_SegCoregNorm2d.m

if prepro.epicoreg.run
    
    fid = fopen('preprolog.txt', 'a');
    fprintf(fid, 'EPI coregistration started at %s\n', datestr(now));
    fprintf('EPI coregistration started at %s\n', datestr(now))
    tic
    
    % Batch Defaults
    % --------------
    % These should be the same as T1 norm step 3
    coregjob{1}.spatial{1}.coreg{1}.estimate.other = {''};
    coregjob{1}.spatial{1}.coreg{1}.estimate.eoptions.cost_fun = 'nmi';
    coregjob{1}.spatial{1}.coreg{1}.estimate.eoptions.sep = [4 2 1]; % default is [4 2]
    coregjob{1}.spatial{1}.coreg{1}.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001]./2; % half of the defaults
    coregjob{1}.spatial{1}.coreg{1}.estimate.eoptions.fwhm = [7 7];
    
    % Coregister to native GM
    gm = spm_select('FPlist', T1DirPath, sprintf('^c1.*%s.*.nii$','t1_mprage'));
    coregjob{1}.spatial{1}.coreg{1}.estimate.ref = {gm};
    
    % The thing to move around, the mean, realigned image
    if ~prepro.epicoreg.runbyrun
        run1Path = fullfile(boldPath, sprintf('%.3d', runDirNums(1)));
        meanfname = dir(fullfile(run1Path, 'meanlf*.nii'));
        if length(meanfname) ~= 1
            error('more than one meanfname')
        end
        coregjob{1}.spatial{1}.coreg{1}.estimate.source = ...
            { fullfile(run1Path, meanfname.name) };
    end
    
    % Even though adding scans for all runs, have to be careful about adding
    %   recursively for all runs at once because there may be bold acquisition
    %   folders that want to leave out
    
    other = {};
    nfiles = 0;
    for ri = 1:length(runDirNums)
        
        % Note that the file names must be the full path and must be a column
        %  vector
        runPath = fullfile(boldPath, sprintf('%.3d', runDirNums(ri)));
        
        % Find all the functional files for this run.
        % Assume that they are prefixed with an 'raf' from realignment step
        scanfnames = dir(fullfile(runPath, 'rlf*.nii'));
        
        % Need the full path, so concatinate the path and the name
        %   Note that fullfile does not add the last "/" so have to add that by
        %   hand
        if prepro.epicoreg.runbyrun
            
            % Get the volume to compare to
            % Will only work if realigned run-by-run as well
            meanfname = dir(fullfile(runPath, 'meanlf*.nii'));
            if length(meanfname) ~= 1
                error('more than one meanfname')
            end
            coregjob{1}.spatial{1}.coreg{1}.estimate.source = ...
                { fullfile(runPath, meanfname.name) };
            
            other = reshape(strcat(runPath, '/', {scanfnames.name}, ',1'), [], 1);
            nfiles = nfiles + length(other);
            
            % Add the scans into the batch
            coregjob{1}.spatial{1}.coreg{1}.estimate.other = other;
            
            % Save the batch
            fname = fullfile(sessPath, sprintf('batchcoregepi_%d.mat', ri));
            save(fname, 'coregjob')
            fprintf(fid, 'Saved coregjob to %s\n', fname);
            
            % Run the batch
            fprintf('\nCoregistering EPIs in %s to structural...', sprintf('%.3d', runDirNums(ri)))
            spm_jobman('run', coregjob)
            
            % Clear just the source & other, all the rest will stay the same
            coregjob{1}.spatial{1}.coreg{1}.estimate.source = {};
            coregjob{1}.spatial{1}.coreg{1}.estimate.other = {};
        else
            other = [other; ...
                reshape(strcat(runPath, '/', {scanfnames.name}, ',1'), [], 1)];
        end
    end
    
    if ~prepro.epinorm.runbyrun
        % Add the scans into the batch
        coregjob{1}.spatial{1}.coreg{1}.estimate.other = other;
        
        % Save the batch
        save(fullfile(sessPath, 'batchcoregepi.mat'), 'coregjob')
        fprintf(fid, 'Saved matlabbatch to %s\n', ...
            fullfile(sessPath, 'batchcoregepi.mat'));
        
        nfiles = length(other);
        
        % Run the batch
        fprintf('\nCoregistering EPIs to structural...')
        spm_jobman('run', coregjob)
        
    else
        fprintf(fid, 'RUN-BY-RUN\n');
        fprintf('RUN-BY-RUN\n');
    end
    
    % Output complete
    fprintf(fid, 'EPI coregistration completed on %d files in %.3g min\n', ...
        nfiles, toc/60);
    fprintf('EPI coregistration completed on %d files in %.3g min\n', ...
        nfiles, toc/60)
    fclose(fid);
    clear coregjob
end

%% EPI Normalization

% Mitchell just writes the normed EPIs from the estimate gotten from the
%   anatomical. Here I think I want to let them normalize based on the run
%   and the scan because of the distortion in the images (that doesn't
%   existi in the T1). I am going to use the parameters that were used on
%   the T1

if prepro.epinorm.run
    
    fid = fopen('preprolog.txt', 'a');
    fprintf(fid, 'EPI Normalization started at %s\n', datestr(now));
    fprintf('EPI Normalization started at %s\n', datestr(now))
    tic
    
    % Note, these parameters should be adjusted for each study! See
    %   FMRI_Prepro_SPM_...docx or the help under the Normalise,
    %   Normalise: Estimate & Reslice button for more details
    if ~prepro.epinorm.runbyrun
        run1Path = fullfile(boldPath, sprintf('%.3d', runDirNums(1)));
        meanfname = dir(fullfile(run1Path, 'meanlf*.nii'));
        if length(meanfname) ~= 1
            error('more than one meanfname')
        end
        meanepiPath = fullfile(run1Path, meanfname.name);
    end
    
    % Estimate options
    % ----------------
    
    % Original options from T1 Norm #4
    estopts.biasreg = 0.0001; % default = 0.0001
    estopts.biasfwhm= 30; % default = 60-75; reduce for monkey brain?
    
    % Rather than normalizing to the McLaren atlas, normalize to the
    %   normalized segmentations for this animal/session (not selecting
    %   bias-corrected ones because the epis aren't bias corrected, but
    %   don't know if it makes a difference).
    gm = spm_select('FPlist', T1DirPath, sprintf('^wc1.*%s.*.nii$','t1_mprage'));
    wm = spm_select('FPlist', T1DirPath, sprintf('^wc2.*%s.*.nii$','t1_mprage'));
    csf = spm_select('FPlist', T1DirPath, sprintf('^wc3.*%s.*.nii$','t1_mprage'));
    estopts.tpm = char({ [gm ',1'], [wm ',1'], [csf ',1'] });
    
    estopts.ngaus = [2 2 2 8]; % default = [2 2 2 4]
    estopts.msk='';
    estopts.regtype='subj';    % turn on affine: subj corresponds to "average sized template" as used by McLaren 2010
    estopts.samp = 2; % default=3; McLaren 2010 used 2
    estopts.warpco = 12; % default = 25; reduce for monkey brain?
    estopts.warpreg = 1; % default =  1

    
    % Write options
    % -------------
    
    % From T1 Norm #6
    defs.write.bb = [-35 -30 -10; 35 55 45];
    defs.write.vox = [0.5 0.5 0.5];
    
    % Even though adding scans for all runs, have to be careful about adding
    %   recursively for all runs at once because there may be bold acquisition
    %   folders that want to leave out
    
    resample = [];
    nfiles = 0;
    for ri = 1:length(runDirNums)
        
        % Note that the file names must be the full path and must be a column
        %  vector
        runPath = fullfile(boldPath, sprintf('%.3d', runDirNums(ri)));
        
        % Find all the functional files for this run.
        % Assume that they are prefixed with an 'raf' from realignment step
        %scanfnames = dir(fullfile(runPath, 'rlf*.nii'));
        scanfnames = spm_select('FPlist', runPath, '^rlf.*\.nii$');
        
        % Need the full path, so concatinate the path and the name
        %   Note that fullfile does not add the last "/" so have to add that by
        %   hand
        if prepro.epinorm.runbyrun
            
            % Get the volume to compare to
            % Will only work if realigned run-by-run as well
            meanfname = dir(fullfile(runPath, 'meanaf*.nii'));
            if length(meanfname) ~= 1
                error('more than one meanfname')
            end
            meanepiPath = fullfile(runPath, meanfname.name);
            
            % Do the Estimation
            % -----------------
            
            % Remember, this is on the coregistered image
            out = spm_preproc(meanepiPath, estopts);
            [sn,isn] = spm_prep2sn(out);
            clear out
            
            % Write forward and inverse deformation files
            % -----
            
            snmatname = [spm_str_manip(meanepiPath,'sd') '_seg_sn.mat'];
            savefields(snmatname, sn);
            savefields([spm_str_manip(meanepiPath,'sd') '_seg_inv_sn.mat'], isn);
            
            % Find the files to normalise
            
            resample = scanfnames;
            nfiles = nfiles + length(resample);
            
            % Write the normalised files
            % --------------------------
            % Inclued the mean EPI!
            resample = [resample; [meanepiPath ',1']];
            spm_write_sn(resample, snmatname, defs.write);
            
            % Clear just the vol & resample, all the rest will stay the same
            clear meanepiPath resample
        else
            % Add this run on to the list
            resample = [resample; scanfnames];
        end
    end
    
    if ~prepro.epinorm.runbyrun
        
        % Do the Estimation
        % -----------------
        
        % Remember, this is on the coregistered image
        out = spm_preproc(meanepiPath, estopts);
        [sn,isn] = spm_prep2sn(out);
        clear out
        
        % Write forward and inverse deformation files
        % -----
        
        snmatname = [spm_str_manip(meanepiPath,'sd') '_seg_sn.mat'];
        savefields(snmatname, sn);
        savefields([spm_str_manip(meanepiPath,'sd') '_seg_inv_sn.mat'], isn);
        
        % Write the normalised files
        % --------------------------
        % Inclued the mean EPI!
        spm_write_sn(char(resample, meanepiPath), snmatname, defs.write);
        
        nfiles = length(resample);
        
    else
        fprintf(fid, 'RUN-BY-RUN\n');
        fprintf('RUN-BY-RUN\n');
    end
    
    % Output complete
    fprintf(fid, 'EPI Normalization completed on %d files in %.3g min\n', ...
        nfiles, toc/60);
    fprintf('EPI Normalization completed on %d files in %.3g min\n', ...
        nfiles, toc/60)
    fclose(fid);
    %clear matlabbatch
end %of normalization

%% Smoothing
% Even if did realignment and normalization run-by-run, can do all runs at
%   once here because smoothing is independent of time

if prepro.smooth.run
    
    fid = fopen('preprolog.txt', 'a');
    fprintf(fid, 'Smoothing started at %s\n', datestr(now));
    fprintf('Smoothing started at %s\n', datestr(now))
    tic
    
    % Parameters
    % ----------
    
    fwhm = 3; % used by Mars et al. and Hutchinson et al.; Mantini used 2mm
    % Cell array of files to restrict smoothing within masks
    %   Referenced from Mitchell_scripts/aamod_smooth.m
    %   % see e.g.
    % https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind1008&L=spm&P=R45120&1=spm&9=A&J=on&K=2&d=No+Match%3BMatch%3BMatches&z=4
    % https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind0901&L=spm&P=R14635&1=spm&9=A&J=on&K=3&d=No+Match%3BMatch%3BMatches&z=4
    masks{1} = fullfile(T1DirPath, 'wGM.nii');
    masks{2} = fullfile(T1DirPath, 'wWM.nii');
    masks{3} = fullfile(T1DirPath, 'wCSF.nii');

    % Get images to reslice
    % ---------------------
    
    % Even though adding scans for all runs, have to be careful about adding
    %   recursively for all runs at once because there may be bold acquisition
    %   folders that want to leave out
    
    subjimgs = [];
    for ri = 1:length(runDirNums)
        
        % Note that the file names must be the full path and must be a column
        %  vector
        runPath = fullfile(boldPath, sprintf('%.3d', runDirNums(ri)));
        
        % Find all the functional files for this run.
        % Assume that they are prefixed with an 'f'
        scanfnames = spm_select('FPlist', runPath, '^wrlf.*\.nii$');
        
        % Need the full path, so concatinate the path and the name
        %   Note that fullfile does not add the last "/" so have to add that by
        %   hand
        subjimgs = [subjimgs; scanfnames];
    end
    
    % Start copying here directly from aamod_smooth.m
    
    % now smooth
    % ----------
    s   = fwhm;
    n   = size(subjimgs,1);
    
    % Smooth the masks, if necessary (files may already exist)
    rmasks{1} = fullfile(T1DirPath, 'rwGM.nii');
    rmasks{2} = fullfile(T1DirPath, 'rwWM.nii');
    rmasks{3} = fullfile(T1DirPath, 'rwCSF.nii');
    if ~isempty(masks)
        fprintf('\nSmoothing masks')
        Vm = cell(1,length(masks));
        Ym = cell(1,length(masks));
        Yms = cell(1,length(masks));
        reslicedmask = cell(1,length(masks));
        for m = 1:length(masks);
            anEPI = deblank(subjimgs(1,:));
            
            [pth, nam] = fileparts(masks{m});
            thismask = spm_select('FPList', pth, nam);
            if ~exist(rmasks{m}, 'file')
                spm_reslice(char(anEPI,thismask),struct('which',1,'mean',0))
            end
            reslicedmask{m} = regexprep(thismask,'(.*)/(.*)','$1/r$2');
            
            Vm{m} = spm_vol(reslicedmask{m});
            Ym{m} = round(spm_read_vols(Vm{m}));
            Yms{m} = nan(size(Ym{m}));
            [BB, vx{m}] = spm_get_bbox(Vm{m});
            spm_smooth(Ym{m}, Yms{m},abs(s./[vx{m}(1) vx{m}(2) vx{m}(3)]));
            Ym{m} = logical(Ym{m});
            fprintf('.')
        end
    end
    
    spm_progress_bar('Init',n,'Smoothing','Volumes Complete');
    fprintf('\nSmoothing EPIs')
    warning off all % don't worry about divide by zero
    for i = 1:n
        Q = deblank(subjimgs(i,:));
        [pth,nam,xt,nm] = spm_fileparts(Q);
        U = fullfile(pth,[sprintf('s%g',s) nam xt nm]); %djm: add smoothing size
        if isempty(masks)
            spm_smooth(Q,U,s);
        else % djm: smooth within masks, but not between them
            VQ = spm_vol(Q);
            YQ = spm_read_vols(VQ);
            YU = YQ;
            for m = 1:length(masks);
                YQs = nan(size(YQ));
                spm_smooth(YQ.*Ym{m},YQs,abs(s./[vx{m}(1) vx{m}(2) vx{m}(3)]));
                YQs = YQs./Yms{m}; % this 'undoes' the zero padding of convolution at mask edges
                YU(Ym{m}) = YQs(Ym{m});
            end
            VQ.fname = U;
            spm_write_vol(VQ,YU);
        end
        spm_progress_bar('Set',i);
        if ~mod(i,50), fprintf('.'); end
    end
    spm_progress_bar('Clear');
    warning on all
    
    
    % --------- Old Batch -------------
    %{
    % Add the scans back into the batch
    matlabbatch{1}.spm.spatial.smooth.data = data;
    
    % Save the batch
    save(fullfile(sessPath, 'batchsmooth.mat'), 'matlabbatch')
    fprintf(fid, 'Saved matlabbatch to %s\n', ...
        fullfile(sessPath, 'batchsmooth.mat'));
    
    % Run the batch
    spm_jobman('run', matlabbatch);
    %}
    
    % Output complete
    fprintf(fid, 'Smoothing completed on %d files in %.3g min\n', ...
        length(subjimgs), toc/60);
    fprintf('Smoothing completed on %d files in %.3g min\n', ...
        length(subjimgs), toc/60)
    fclose(fid);
    %clear matlabbatch
end



%% Finishing

% Set the path back to the original path before finishing
%Shouldn't really need this for CCV
%cd(origWorkPath)
fprintf('**********************\nDone preprocessing %s\n**********************\n', ...
    sessPath)