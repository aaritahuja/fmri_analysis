%<subjnum> input as a number

% Code taken from SPM mailing list
% --- from here:
% https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=SPM;e4a3a9e2.1207 ---%

function realignment_motion_check(subjnum, subjPath)

%subjPath = fullfile('/gpfs/data/tdesroch/pSEQ/subjects', int2str(subjnum));
boldPath = fullfile(subjPath, 'bold');
rundirs = {'001'}; % If realignment was done all at once, the information will be in a single file located in the folder for the first run.
%If done one at a time, each folder will have its own file, in which case
%see below

%Added loop for runbyrun graphs
%rundirs = {'001' '002' '003' '004' '005' '006' '007' '008' '009' '010' '011' '012' '013' '014'};

for runi = 1:length(rundirs) % needed for RUN by RUN realignment
    
    cd(fullfile(boldPath, sprintf('%.3d',runi)))

    motion_file = dir('*rp*');
    motion_parameters = load(motion_file.name);

    scaleme = [-3 3];

    printfig = figure;
    set(printfig, 'Name', ['Motion parameters: subject ' num2str(subjnum) ], 'Visible', 'on');
    subplot(2,1,1);
    plot(motion_parameters(:,1:3));
    grid on;
    legend('x','y','z');
    ylim(scaleme);  % enable to always scale between fixed values as 
                     % set above
    title([ sprintf('Data from %s\n\n' , pwd) 'Translation (in mm)']);
           % 'interpreter', 'none');
    subplot(2,1,2);
    plot(motion_parameters(:,4:6)*180/pi);
    grid on;
    legend('x','y','z');
    ylim(scaleme);   % enable to always scale between fixed values as 
                      % set above
    title('Rotations (in dg)', 'interpreter', 'none');
    mydate = date;
    motname = [pwd filesep 'motion_sub_' sprintf('%s', num2str(subjnum)) '_' mydate '.png'];


end