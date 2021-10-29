function td_biasCorrect(T1Path)

% <T1Path> is the full string path name of the file to de-bias. Will output
%   a file prepended with 'm'

% Taking a section of td_prepro 5/18/17

% Last modified 05/18/17 by Theresa M. Desrochers

%% Paths

if ismac
    atlasPath =  '/Users/Aaru/Documents/MRI_monkey/atlasfiles';
else
    atlasPath = 'C:/Users/lab/Documents/MRI_monkey/atlasfiles';
end
templatePath = fullfile(atlasPath, '112RM-SL_T1.nii');

[T1DirPath, fname, ext] = fileparts(T1Path);
% fname doesn't have extension on it
T1name = [fname ext]; 
mT1name = sprintf('m%s',T1name);
mT1Path = fullfile(T1DirPath, mT1name);


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
    fprintf('\nWriting bias-corrected T1...')
    spm_preproc_write(sn, writeopts);
else
    fprintf('\nBias-corrected T1 found.')
end