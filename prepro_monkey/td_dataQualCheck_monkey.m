function final_table = td_dataQualCheck_monkey(sessPath, runDirs, fprefix, varargin)
% Merge of td_dataQualCheck original version made by TMD and the modified
% version by AA (2019).
% This version will now save all the figures separately, as well as output
% the table. 
% Inputs:
% <boldPath> full path string of the upper level directory
% <runDirs> cell array of run folders
% <fprefix> string file prefix, usually f, af, raf, or wraf

% Based on rawQuant.m by CG
% Last modified 12/22/2019 by NYR
tic

figurePath = fullfile(sessPath, 'qualcheck');

acquisitions = strings([length(runDirs),1]);
variance = strings([length(runDirs),1]);
drifts = strings([length(runDirs),1]);
max_variance = strings([length(runDirs),1]);
stdev = strings([length(runDirs),1]);
default_threshold = strings([length(runDirs),1]);
movement = strings([length(runDirs),1]);
boldPath = fullfile(sessPath, 'bold')
SNR = td_CheckSNR_monkey(boldPath, runDirs);

for ri = 1:length(runDirs)
    %fidx = ~cellfun(@isempty, strfind(subjFiles, sprintf('_run%d', runnums(ri))));
    
    rundir = fullfile(boldPath, runDirs{ri} );
    
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
    set(gcf,'visible','off')
    %tsdiffana
    % default position [0.3256    0.0084    0.3750    0.9105]
    hF(1) = figure('Units', 'normalized', 'Position', [0,0,.25,.9]);
    tsdiffana_sub(dirnames, 0, hF(1));
    %Extract values
    fig1 = gcf;
    axObjs = fig1.Children;
    axchild = allchild(axObjs); 
    %First figure: Outputs min and max values
    %AA  had done a lot to find these vals, but CCV hates opengl
    %This makes the value sbe different so changed them to match in CCV
    if isunix
        ydata = axchild(5);
    else 
        ydata = axchild(5);
    end
    %ydata = axchild(14);
    y_conv = ydata{1};
    y_vals = y_conv.YData;
    first_graph_max = max(y_vals);
    first_graph_min = min(y_vals);
    variance_range = [num2str(first_graph_min), ' : ', num2str(first_graph_max)];
    variance(ri, 1) = variance_range;
    
    %Third figure: Check for drift
    %Same problem here
    if isunix
        ydata = axchild(3);
    else
        ydata = axchild(3);
    end
    %ydata = axchild(12);
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
    if isunix
        ydata = axchild(2);
    else 
        ydata = axchild(2);
    end
    %ydata = axchild(11);
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
    if isunix
        sdev_graph = axObjs(26);
    else
        sdev_graph = axObjs(26);
    end
    %sdev_graph = axObjs(35);
    %These std_YData values below are the ones we want to save in a .text
    %file for each run to make the regressor
    std_YData = sdev_graph.Children(end).YData';
    %Save to a .txt file for this variable with the sessnum attached
    fileparts = regexp(regexp(boldPath,filesep,'split'), '\d+(\.)?(\d+)?','match');
    sessnum = sprintf('0%d',str2double([fileparts{:}]));
    std_filename = strcat('std_',sessnum);
    file = [rundir filesep std_filename];
    save(file, 'std_YData', '-ASCII')
    
    %Alignment axis range
    if isunix
        alignment_axis = axObjs(25);
    else
        alignment_axis = axObjs(25);
    end
    %alignment_axis = axObjs(34);
    max_motion = max(alignment_axis.YLim);
    min_motion = min(alignment_axis.YLim);
    motion_range = [num2str(max_motion), ' : ', num2str(min_motion)];
    movement(ri, 1) = motion_range;
    %Save figure
    saveas(gcf, strcat(figurePath, '/run', num2str(ri), '_figure2.png'))
    
    % art_movie
    %hF(3) = figure('units', 'normalized', 'position', [.65, 0, .35, .8]);
    %art_movie_sub(1,dirnames,hF(3));
    %This saves the figure, however because this is in fact a movie this is
    %not very useful. If the data is too bad, it might be best to review
    %the movie using the original version of this code.
    %saveas(gcf, strcat(figurePath, '/run', num2str(ri), '_figure3.png'))
    %% Commented out the parts of the code that actually oputput and cycle through figs.
    % Make an interact flag so that can interact if want.
%     fprintf('%s Run %d: %d files\n', boldPath, runDirs{ri}, length(fnames))
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
%% Want the table with info as output
final_table = table(runDirs', acquisitions, SNR, variance, drifts, max_variance, stdev, default_threshold, movement);
toc
   