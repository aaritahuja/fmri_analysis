
%<subjnum> input as a number

% Code taken from SPM mailing list
% --- from here:
% https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=SPM;e4a3a9e2.1207 ---%

function db_realignment_motion_check(subjnum)

subjPath = fullfile('/gpfs/data/tdesroch/pSEQ/subjects', int2str(subjnum));
boldPath = fullfile(subjPath, 'bold');
rundirs = {'001'}; % Checking for first run
runbyrun=0;

if runbyrun ==0
 rundirs = {'001'};
else
%Added loop for runbyrun graphs
rundirs = {'001' '002' '003' '004' '005' '006'};
end


for runi = 1:length(rundirs) % needed for RUN by RUN realignment
    
    cd(fullfile(boldPath, sprintf('%.3d',runi)))
    
    motion_file = dir('*rp*');
    motion_parameters = load(motion_file.name);
    
    
    
    scaleme = [-4 4];
    
    printfig = figure;
    set(printfig, 'Name', ['Motion parameters: subject ' num2str(subjnum) ], 'Visible', 'on');
    subplot(2,1,1);
    h = plot(motion_parameters(:,1:3));
    grid on;
    ylim(scaleme);  % enable to always scale between fixed values as
    % set above
    
    hold on
    line(xlim, [3 3],'LineStyle','--','LineWidth',1, 'Color', 'k' )
    line(xlim, [-3 -3],'LineStyle','--','LineWidth',1, 'Color', 'k' )
    legend(h,'x','y','z','Location','southwest');
    
    title([ sprintf('Data from %s\n\n' , pwd) 'Translation (in mm)']);
    % 'interpreter', 'none');
    subplot(2,1,2);
    plot(motion_parameters(:,4:6)*180/pi);
    grid on;
    legend('x','y','z','Location','southwest');
    ylim(scaleme);   % enable to always scale between fixed values as
    % set above
    title('Rotations (in dg)', 'interpreter', 'none');
    
    
    mydate = date;
    motname = [pwd filesep 'motion_sub_' sprintf('%s', num2str(subjnum)) '_' mydate '.png'];
    saveas(gcf,motname)
end









%
% for runi = 1:length(rundirs) % needed for RUN by RUN realignment
% %cd(fullfile(boldPath, '001')) % needed for realign normal
%
%
% %this is the dir where rp file is located- contains data for all 6 runs
%
% rp = spm_load(spm_select('List','rp*.txt'));
%
% end