
if ismac
    % To differentiate between different mac computers, use the following:
    % [idum,hostname]= system('hostname');
    addpath('~/Documents/MATLAB/spm12')
elseif ispc
    addpath('~/Documents/MATLAB/spm12')
elseif isunix
    addpath('/gpfs/runtime/opt/spm/spm12')
else
    error('Need to add path to spm12 here')
end