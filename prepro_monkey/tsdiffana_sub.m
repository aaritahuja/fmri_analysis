function [td, globals, slicediff, imgs] = tsdiffana_sub(imgs, vf, fg, pflags)
% wrapper and plotter for timediff function
%
% imgs    - string list of images
% vf      - non zero if difference images required
% fg      - figure to plot results (spm figure default)
% pflags  - plot flags - 0 or more of 'r'  - plot realignment params

% tsdiffana_sub.m by CG & TMD to not create it's own figure Nov 2012.

if nargin < 1
    imgs = [];
end
if isempty(imgs)
    imgs = spm_select;
end
if nargin < 2
    vf = 0;
end
noplot = 0;
if nargin < 3
    fg = [];
elseif isempty(fg)
    noplot = 1;
end
if nargin < 4
    pflags = '';
end

[p f e] = fileparts(imgs(1,:));
if vf
    flags = 'mv';
else
    flags = 'm';
end

%keyboard

[td globals slicediff] = timediff(imgs,flags);
tdfn = fullfile(p,'timediff.mat');
save(tdfn, 'imgs', 'td', 'globals', 'slicediff');


if ~noplot
    tsdiffplot(tdfn,fg, pflags);
end
return

