function td_realign(anatPath, imgPaths)

% <anatPath> string of anatomy folder that imgPaths can be found in
% <imgPaths> cell array of full paths of images that want to re-aling to
%   each other. Path should start as if adding on to the anatPath folder. 
%   Alignment will be done to first image.

% Taking a section of td_prepro 5/18/17

% Last modified 05/18/17 by Theresa M. Desrochers

origWorkPath = cd;
cd(anatPath)

%% Realignment

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
%   Load an example image to get the dimensions - this is only done for
%   first image!
v = spm_vol(fullfile(anatPath, imgPaths{1}));
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
%   1 = Register to mean: A two pass procedure is used in order to register
%   the images to the mean of the images after the first realignment
%   0 = Register to first image: as will weight the images for the average
%   later, do this now.
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;

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


data = reshape(strcat(fullfile(anatPath, imgPaths), ',1'), [], 1);

% Add the scans into the batch
matlabbatch{1}.spm.spatial.realign.estwrite.data = {data}';

% Save the batch
save(fullfile(anatPath, 'batchT1realalign.mat'), 'matlabbatch')

nfiles = length(data);

% Run the batch
spm_jobman('run', matlabbatch);

% Set the path back to the original path before finishing
cd(origWorkPath)

% Output complete
fprintf('Realignment completed on %d files in %.3g min\n', ...
    nfiles, toc/60)
clear matlabbatch
