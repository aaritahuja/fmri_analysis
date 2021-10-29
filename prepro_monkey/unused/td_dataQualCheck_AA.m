function final_table = td_dataQualCheck_AA(MRI_path, subjnum, runnums, fprefix)

% Need to run td_file_sort first!
% Inputs:
% <subjnum> subject number
% <runnums> array of run numbers
% <fprefix> string file prefix, usually f, af, raf, or wraf
% Example: td_dataQualCheck(102, 1:6, 'f') 

% Based on rawQuant.m by CG
% Last modified 11/10/15 TMD
% Example: td_dataQualCheck(102, 1:6, 'f') 

subjPath = fullfile(MRI_path, sprintf('%d', subjnum));
% subjFiles = dir(subjPath);
% subjFiles = {subjFiles.name};
boldPath = fullfile(subjPath, 'bold');
figurePath = fullfile(subjPath, 'preprocessing figures');


acquisitions = strings([length(runnums),1]);
variance = strings([length(runnums),1]);
drifts = strings([length(runnums),1]);
max_variance = strings([length(runnums),1]);
stdev = strings([length(runnums),1]);
default_threshold = strings([length(runnums),1]);
movement = strings([length(runnums),1]);

for ri = 1:length(runnums)
    %fidx = ~cellfun(@isempty, strfind(subjFiles, sprintf('_run%d', runnums(ri))));
    
    rundir = fullfile(boldPath, sprintf('%.3d',runnums(ri)) );
    
    % Get the functional files for the current bold directory
    % May want to fix to pre-pend the files with f
    dirdata = dir(fullfile(rundir, [fprefix '*.nii']));
    n_acq = num2str(length(dirdata));
    acquisitions(ri, 1) = n_acq; 
    
    % Convert struct to character array
    fnames = char(dirdata.name);
    % Make it the full file name
    dirnames = strcat(rundir, '/', fnames);
    
    % Three figures
    
    % tsdiffana
    % default position [0.3256    0.0084    0.3750    0.9105]
    hF(1) = figure('Units', 'normalized', 'Position', [0,0,.25,.9]);
    tsdiffana_sub(dirnames, 0, hF(1));
    %Extract values
    fig1 = gcf;
    axObjs = fig1.Children;
    axchild = allchild(axObjs); 
    %First figure: Outputs min and max values
    ydata = axchild(14);
    y_conv = ydata{1};
    y_vals = y_conv.YData;
    first_graph_max = max(y_vals);
    first_graph_min = min(y_vals);
    variance_range = [num2str(first_graph_min), ' : ', num2str(first_graph_max)];
    variance(ri, 1) = variance_range;
    
    %Third figure: Check for drift
    ydata = axchild(12);
    ydata = ydata{1};
    ydata = ydata.YData;
    xlist = 1:length(dirdata);
    b = corrcoef(xlist, ydata);
    coefficient = b(1,2);
    if coefficient <= -0.3
        drift = 'Negative';
    elseif (coefficient >= -0.3) && (coefficient <= -0.2)
        drift = 'Slight negative';
    elseif (coefficient >= -0.2) && (coefficient <= 0.2)
        drift = 'None';
    elseif (coefficient >= 0.2) && (coefficient <= 0.3)
        drift = 'Slight positive';
    elseif coefficient >= 0.3
        drift = 'Positive';
    end
    drifts(ri, 1) = drift;
        
    %Last figure: Outputs max value
    ydata = axchild(11);
    y_conv = ydata{1};
    y_vals = y_conv.YData;
    last_graph_max = num2str(max(y_vals));
    max_variance(ri, 1) = last_graph_max;
    %Save figure
    saveas(gcf, strcat(figurePath, '/run', num2str(ri), '_figure1.png'))
    
    % art_global
    hF(2) = figure('units', 'normalized', 'position', [.25, 0, .4, .8]);
    art_global_sub(dirnames, hF(2));
    %Extract values
    fig2 = gcf;
    axObjs = fig2.Children;
    object_handles = findobj(fig2, 'type', 'uicontrol');
    %Current threshold
    threshold_table = object_handles(18);
    threshold = threshold_table.String;
    default_threshold(ri, 1) = threshold;
    %StdDev
    sdev_table = object_handles(20);
    sdev = sdev_table.String;
    stdev(ri, 1) = sdev;
    %Alignment axis range
    alignment_axis = axObjs(34);
    max_motion = max(alignment_axis.YLim);
    min_motion = min(alignment_axis.YLim);
    motion_range = [num2str(max_motion), ' : ', num2str(min_motion)];
    movement(ri, 1) = motion_range;
    %Save figure
    saveas(gcf, strcat(figurePath, '/run', num2str(ri), '_figure2.png'))
    
    
    % art_movie
    hF(3) = figure('units', 'normalized', 'position', [.65, 0, .35, .8]);
    art_movie_sub(1,dirnames,hF(3));
    saveas(gcf, strcat(figurePath, '/run', num2str(ri), '_figure3.png'))
    
%     fprintf('Subj %d Run %d: %d files\n', subjnum, runnums(ri), length(fnames))
%     inp = input('Next run & close figs (n), next run & open new figs (o), exit (x)? ', 's');
%     switch inp
%         case 'n'
%             close(hF);
%         case 'o'
%             clear hF
%         case 'x'
%             break
%         otherwise
%             error('unrecognized option')
%     end
    close(hF)
end
final_table = table(acquisitions, variance, drifts, max_variance, stdev, default_threshold, movement);
   