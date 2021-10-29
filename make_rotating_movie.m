function make_rotating_movie(subj)

ROI_file = fullfile('/Users/Aaru/Documents/MRI_data', num2str(2), 'models/localizer/higher_order_motion.mat');
ROI = load(ROI_file);

ROI_X = ROI.higher_order_motion_coordinates(:, 1);
ROI_Y = ROI.higher_order_motion_coordinates(:, 2);
ROI_Z = ROI.higher_order_motion_coordinates(:, 3);

cd(fullfile('/Users/Aaru/Documents/MRI_data', num2str(subj), 'models/localizer'))
load('SPM.mat');
t = zeros(length(ROI_X), 1);

for i = 1:length(ROI_X)
    coords = [ROI_X(i); ROI_Y(i); ROI_Z(i)];
    [currbeta, currSE, currCI] = extractSPMData(SPM, coords, contrast);
    curr_t = currbeta/currSE;
    t(i, 1) = curr_t;
end

anatomical_epi_l = niftiread(fullfile('/Users/Aaru/Documents/MRI_data', num2str(subj), 'anatomy/wsurface-volume-left.nii'));
anatomical_l = anatomical_epi_l >= 1;
[xl, yl, zl] = ind2sub(size(anatomical_epi_l), find(anatomical_l == 1));
cl = anatomical_epi_l(anatomical_l);

anatomical_epi_r = niftiread(fullfile('/Users/Aaru/Documents/MRI_data', num2str(subj), 'anatomy/wsurface-volume-right.nii'));
anatomical_r = anatomical_epi_r >= 1;
[xr, yr, zr] = ind2sub(size(anatomical_epi_r), find(anatomical_r == 1));
cr = anatomical_epi_r(anatomical_r);

x = vertcat(xl, xr);
y = vertcat(yl, yr);
z = vertcat(zl, zr);
c = vertcat(cl, cr);

figure
hold on
set(gca,'visible','off')
set(gcf, 'Position',  [100, 100, 1200, 900])

view(54, 50)
roi = scatter3(ROI_X, ROI_Y, ROI_Z, 40, t, 'filled', 's', 'MarkerFaceAlpha', 0.2, 'MarkerFaceColor', [0.65 0.65 0.65]);
colormap hot
colormax = max(t);
colormin = 0;
caxis([colormin colormax])
brain = scatter3(x, y, z, 40, 'filled', 's', 'MarkerFaceAlpha', 1, 'MarkerFaceColor', [0.5 0.5 0.5]);

axis vis3d
for i = 1:361
camorbit(1,0,'data',[0 0 1])
F(i) = getframe(gcf);
drawnow
end

% create the video writer with 1 fps
writerObj = VideoWriter('with_roi.mp4', 'MPEG-4');
writerObj.FrameRate = 60;

% open the video writer
open(writerObj);
% write the frames to the video
for a=1:length(F)
    % convert the image to a frame
    frame = F(a) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);