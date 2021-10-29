function td_write_sn(niiPath, snPath)

% td_write_sn(niiPath, snPath)
% <niiPath> the full string path of the file to be transformed
% <snPath> the full string path of the sn.mat file to use to transform the
%   file.
% New file will be saved with a 'w' prepended to the file name in the
%   original directory.
% Need to have SPM in your path, uses an SPM function

% 3/8/17 TMD

% Default write options, taken from td_prepro
defs.write.bb = [-35 -30 -10; 35 55 45];
defs.write.vox = [0.5 0.5 0.5];

spm_write_sn(niiPath, snPath, defs.write);