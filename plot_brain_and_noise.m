%% Makes a 3D plot of a functional volume of the brain as well as the regions that would be considered noise

function [brain_width, img_contrast] = plot_brain_and_noise (MRI_path, subject, run, varargin)

%Get paths and filenames - Different for human (int) and monkey (char)
if isinteger(subject)
    subjectBoldPath = fullfile(MRI_path, num2str(subject), 'bold');
elseif ischar(subject)
    if  nargin > 0
        session = varargin{1};
    else
        error('Session number not specified')
    end
    subjectBoldPath = fullfile(MRI_path, subject, session, 'bold');
end

runPath = fullfile(subjectBoldPath, run);
if isinteger(subject)
    files = dir(fullfile(runPath,'f*.nii'));
elseif ischar(subject)
        files = dir(fullfile(runPath,'lf*.nii'));

end
randomidx = randi(length(files));
filename = fullfile(runPath, files(randomidx).name);

%Get nifti file
[epi,pixdim,rotate,dtype] = readnifti(filename);

%Threshold out brain as being above a certain (arbitrary) voxel value
if isinteger(subject)
    brain = epi >= 400;
elseif ischar(subject)
    brain = epi >= 220;
end

%convert to xyz coordinates
[x y z] = ind2sub(size(epi), find(brain == 1));
brain_width = range(x);

collapsed_epi = permute(epi, [2 1 3]);
img_contrast = range(collapsed_epi(:));

%use original values as colorscheme
color = epi(brain);

%initiate noise matrix of same size
noise = epi < 0;

%set noise ranges
noise(2:6,2:6,2:6) = 1;
noise(end-6:end-1,end-6:end-1,end-6:end-1) = 1;
noise(end-6:end-1,end-6:end-1,2:6) = 1;
noise(end-6:end-1,2:6,end-6:end-1) = 1;
noise(2:6,end-6:end-1,end-6:end-1) = 1;
noise(2:6,2:6,end-6:end-1) = 1;
noise(end-6:end-1,2:6,2:6) = 1;
noise(2:6,end-6:end-1,2:6) = 1;

%convert to xyz coordinates
[nx ny nz] = ind2sub(size(noise), find(noise == 1));

max_color = max(max(max(epi)));
min_color = min(min(min(epi)));
%plot
draw = figure;
scatter3(x, y, z, 50, color, 'filled', 'square');
colormap gray
caxis([min_color max_color]);
colorbar;
hold on
scatter3(nx, ny, nz, 50);
view(126, 38)
view(-0.6, 90)
figure(draw)
view(220, 36)


