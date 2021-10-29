function PSTH  = getpsth(xY,SPMfile,window,bin,events,sessions,adjcurr,BETAHACK)
% GETPSTH Extract data from specified events and sessions using FIR model
%   PSTH = GETPSTH(XY,SPMFILE,WINDOW,BIN,EVENTS,SESSIONS,ADJCURR,BETAHACK)
%   returns a structure containing data for plotting a peristimulus time
%   histogram. This differs from the machinery in SPM in that it can
%   extract data for several events at once and across sessions.
%
%   SPM MUST BE IN THE MATLAB PATH. SHOULD WORK WITH BOTH SPM2 AND SPM5 BUT
%   IT HAS NOT BEEN FORMALLY TESTED WITH SPM5.
%
%   Assumes the same number of columns per session in the design matrix.
%   Inputs: xY        Can be either a structure similar to that from
%                     spm_regions, or an [x y z] vector specifying a point.
%           SPMfile   Name of SPM file.
%           window    Length of psth to model.
%           bin       Size of bin for psth in secs (Usually the same as the
%                     TR.
%           events    Event numbers to extract. This corresponds to the
%                     order listed in SPM.Sess.U.
%           sessions  Session numbers to extract
%           adjcurr   Flag of whether to adjust for other event-types in
%                     the same session.
%           BETAHACK  Value to limit betas in case of poor FIR estimates.
%
%   NOTE NOTE NOTE: Cannot deal with sessions that have different event
%   types. It assumes that all event types are present in all sessions!!!
%
%   The function can be called non-interactively by specifying all the
%   inputs.
%
%   Example
%
%      PSTH = getpsth([-21 -69 57],fullfile(pwd,'SPM.mat'),32,2.1,2,1,1)

%   This would extract data from the SPM.mat file in the current directory,
%   at the specfied point, with a window of 32 secs, bin of 2.1 secs, for
%   event 2, session 1, and full adjustments as in spm_graph.
%
%   Output: PSTH      Structure containg the extracted data.
%           SPMfile   same as input.
%           xY        Extracted VOI info
%           window    Specified window length in secs
%           bin       Specified psth bin. Note that the actual window and
%                     bin size may differ from the specified size and is
%                     based on the the requirements of setting up the FIR
%                     basis.
%           adjcurr   Flag of whether to adjust for other columns in the
%                     current session. (i.e., columns other than the event-
%                     type specified.
%           pst       Peristimulus time vector
%           events    Structure of data for each event type.
%             number    Number of the extracted event, based on SPM.Sess.U
%             name      Name of the extracted event.
%             sessions  Extracted sessions for each event. Same as input.
%                       order corresponds to columns in psth.
%             psth      The extracted psth for each session in columns
%             sem       Standard error of the mean for each time point
%             pci       Confidence interval for each time point.
%
%   Example
%      To plot the PSTH with 95% confidence interval errorbars
%
%         errorbar(PSTH.pst,PSTH.events(1).psth,PSTH.events(1).pci)
%
%      To plot the first 2 events, with S.E. errorbars, and a legend.
%      Use n to index the structure so to avoid rewriting the command each
%      time.
%
%         n = 1;
%         errorbar(PSTH.pst,PSTH.events(n).psth,PSTH.events(n).sem)
%         hold on
%         n = 2;
%         errorbar(PSTH.pst,PSTH.events(n).psth,PSTH.events(n).sem)
%         legend(PSTH.events(1).name,[PSTH.events(1).name, ' S.E.'],...
%                PSTH.events(2).name,[PSTH.events(2).name, ' S.E.'])

%      Note that each plotted line gets 2 entries in the legend because matlab
%      insists on having a legend entry for both a line and its errorbars.
%
%      Multiple sessions can be averaged. PSTH.events(n).psth is organized
%      as rows = time bins and cols = sessions. In order to average across
%      sessions, and not time bins, the matrix must be transposed, then
%      averaged, then transposed back. To plot the first event for the
%      average psth across sessions:
%
%         plot(PSTH.pst,mean(PSTH.events(1).psth')')

%   Based on spm_graph and the extract_psth code from Alexa Morcom.

%  Author: Darren Gitelman
%  $Id: getpsth.m,v 1.23 2008-04-11 14:31:59-05 drg Exp drg $

% Check for spm defaults
% --------------------------------------------------------------
var = whos('global');
var = str2mat(var.name);
if ~strmatch('defaults',var)
    error('SPM defaults not present. Please start spm or run spm_defaults first.');
end

% Find or create the interactive window
% --------------------------------------------------------------
Finter = spm_figure('Findwin','Interactive');
if isempty(Finter)
    Finter = spm('CreateIntWin');
end
spm_input('!SetNextPos',1);

% PSTH structure
% --------------------------------------------------------------
PSTH = struct(...
    'SPMfile',  '',...
    'adjcurr',  [],...
    'bin',      [],...
    'pst',      [],...
    'window',   [],...
    'xY',       [],...
    'events',   struct(...
    'name',     '',...
    'number',   [],...
    'pci',      [],...
    'psth',     [],...
    'sem',      [],...
    'sessions', [],...
    'bad',      struct(...
    'pci',      [],...
    'psth',     [],...
    'sem',      [],...
    'sessions', []...
    )...     % close "bad" struct
    )...     % close "events" struct.
    );     % close "PSTH" struct

% CHECK INPUT ARGUMENTS
% --------------------------------------------------------------
if nargin >= 2
    [pth, fn, ext] = fileparts(SPMfile);
    cd(pth)
    load('SPM.mat')
else
    SPM = [];
end

if nargin < 1
    xY = [];
    [PSTH.xY, SPM, SPMfile] = my_regions(xY,SPM);
else
    if ~isstruct(xY)
        if numel(xY) == 3;
            xY.xyz = xY;
        end
    end
    [PSTH.xY, SPM, SPMfileTmp] = my_regions(xY,SPM);
end

if ~exist('SPMfile','var')
    SPMfile = [];
end
if isempty(SPMfile)
    SPMfile = SPMfileTmp;
end;

if isempty(SPMfile)
    switch spm('ver')
        case 'SPM2'
            SPMfile = spm_get(1,'SPM.mat','Select SPM.mat');
        case 'SPM5'
            SPMfile = spm_select(1,'^SPM\.mat$','Select SPM.mat');
        otherwise
            error('Unknown version of SPM.');
    end
    PSTH.SPMfile = SPMfile;
else
    PSTH.SPMfile = SPMfile;
end

[pth fn ex]=fileparts(SPMfile);
if ~isempty(pth)
    cd(pth)
else
    return
end
load(SPMfile);
spm_input('!SetNextPos',1);

if nargin < 4
    % do not check if nargin < 3 so we'll get window and bin together.
    % -----------------------------------------------------------------
    PSTH.window = spm_input('length of psth in seconds','!+1','e',32);
    PSTH.bin    = spm_input('size of each bin in seconds','!+1','e',SPM.xY.RT);
else
    PSTH.window = window;
    PSTH.bin    = bin;
end

if nargin < 5

    % the code assumes that all sessions have at least some representation
    % of all event types so we use Sess(1) as an example
    % --------------------------------------------------------------------
    str = [];
    k = 0;
    str = sprintf('The event-types are:');

    % NOTE: THIS DOES NOT DEAL WITH PARAMETRIC EFFECTS!!!  this may be a
    % silly statment though, since it's not clear that one can plot
    % parameteric effects.
    % MORE IMPORTANTLY IT WILL NOT WORK IF THE DIFFERENT SESSIONS DO NOT
    % ALL HAVE THE SAME EVENT TYPES!!!

    for i = 1:size(SPM.Sess(1).U,2)
        str = str2mat(str,sprintf('[%d]  %s',i,char([SPM.Sess(1).U(i).name{1}])));
    end

    % show the event list as a help dialog
    % --------------------------------------------------------------
    delete(findobj('Tag','eventWindow'));
    h = helpdlg(str,'Event-Types');
    set(h,'Tag','eventWindow');
    figure(Finter);
    events = spm_input('enter event-types','!0','e',[]);
    try
        delete(h)
    catch
    end
end
if nargin < 6
    sessions = spm_input(sprintf('enter sessions [1..%d]',size(SPM.Sess,2)),...
        '!+1','e',[]);
end

if nargin < 7
    % Adjust for in session effects
    % Yes should be the default.
    % -----------------------------------------------------------------
    PSTH.adjcurr = spm_input('Adjust for other in-session effects?','!+1','y/n',[1 0],1);
else
    PSTH.adjcurr = adjcurr;
end


% This option is hidden from users using the GUI. Advanced users will deal
% with this by changing value or using command line.
% --------------------------------------------------------------------
if nargin == 7
    % BETAHACK This variable sets a limit to the beta value. I found that
    % occasionally the least squares equations return an incorrect or
    % wildly high value if there is no effect at a voxel or if only a few
    % events are being used to estimate the effect (2 or less). This value
    % tries to put these beta values in a BAD value field. For my data 100
    % seemed to work. Your data may be different. For no restrictions, just
    % set to Inf.
    % -----------------------------------------------------------------
    PSTH.BETAHACK = BETAHACK;
else
    PSTH.BETAHACK = 100;
end
%--------------------------------------------------------
Now do the work
%--------------------------------------------------------
for n = 1:length(events)
    ev = events(n);
    PSTH.events(n).number = ev;

    % The function assumes that all event-types are
    % present in all sessions.
    % -----------------------------------------------------------------
    PSTH.events(n).name     = SPM.Sess(1).U(ev).name;

    for m = 1:length(sessions)
        ses = sessions(m);

        xBF          = SPM.xBF;
        U            = SPM.Sess(ses).U(ev);

        % event vector in microtime. The event vector is always first in
        % U.u and can be followed by vectors for parametric effects.
        % -----------------------------------------------------------
        U.u          = U.u(:,1);
        xBF.name     = 'Finite Impulse Response';
        xBF.order    = round(PSTH.window/PSTH.bin);
        xBF.length   = xBF.order*PSTH.bin;
        xBF          = spm_get_bf(xBF);
        xBF.bin      = xBF.length/xBF.order;

        if isempty(PSTH.pst)
            j = round(xBF.length/xBF.bin);
            PSTH.pst = [1:j] * xBF.bin - xBF.bin/2;
        end
        X            = spm_Volterra(U,xBF.bf,1);
        k            = SPM.nscan(ses);
        % FIR basis in macrotime
        % -------------------------------------------------------------
        X           = X((0:(k - 1)) * SPM.xBF.T + SPM.xBF.T0 + 32,:);

        % jX will index the appropriate rows of the design. It will
        % increment for different sessions for example.
        % -------------------------------------------------------------
        jX          = SPM.Sess(ses).row;


        if PSTH.adjcurr

            % if adjcurr is true that means the user wants to adjust for other
            % columns in the same session as our chosen event. In that case iX
            % will just be the columns in the original design associated with
            % the chosen event. These are then set to 0 but all the other
            % columns for the session will be present and hence adjusted
            % for
            % ---------------------------------------------------------
            iX = SPM.Sess(ses).col(SPM.Sess(ses).Fc(ev).i);
        else

            % if adjcurr is false then we don't want to adjust for other columns in
            % the same session. So then iX is now equal to all the columns
            % in a session which are then set to 0.
            % ---------------------------------------------------------
            iX = SPM.Sess(ses).col([SPM.Sess(ses).Fc.i]);
        end

        % iX0 represents all the columns of the design but with the columns
        % of the selected event type removed, i.e., remove iX from iX0
        % -------------------------------------------------------------
        iX0         = 1:size(SPM.xX.X,2);
        iX0(iX)     = [];

        % get the selected rows and columns of the design. Note that this
        % moves the effects of interest to the head of the class (i.e., the
        % first set of columns of the matrix.
        % -------------------------------------------------------------
        X           = [X SPM.xX.X(jX,iX0)];

        % whiten the design
        % -------------------------------------------------------------
        X           = SPM.xX.W(jX,jX)*X;

        % add the filter matrix
        % -------------------------------------------------------------
        X           = [X SPM.xX.K(ses).X0];

        % Re-estimate to get PSTH and CI
        %------------------------------------------------------
        CI          = 1.6449;     % = spm_invNcdf(1-0.05);
        dt          = U.dt;
        j           = xBF.order;
        xX          = spm_sp('Set',X);
        pX          = spm_sp('x-',xX);
        betas       = pX*PSTH.xY.u(jX);
        res         = spm_sp('r',xX,PSTH.xY.u(jX));
        df          = size(X,1) - size(X,2);
        bcov        = pX*pX'*sum(res.^2)/df;
        betas       = betas(1:j)/dt;
        sem         = sqrt(diag(bcov(1:j,(1:j))))/dt;
        pci         = CI*sem;

        % this is the hack to trap for wildly incorrect betas. This seems
        % to happen if a voxel has no effects or there are very few events
        % The cutoff chosen is liberal for my data but may not be correct
        % for yours
        % -------------------------------------------------------------
        if max(abs(betas)) < PSTH.BETAHACK
            PSTH.events(n).pci      = [PSTH.events(n).pci,      pci];
            PSTH.events(n).psth     = [PSTH.events(n).psth,     betas];
            PSTH.events(n).sem      = [PSTH.events(n).sem,      sem];
            PSTH.events(n).sessions = [PSTH.events(n).sessions, sessions(m)];
        else
            PSTH.events(n).bad.pci       = [PSTH.events(n).bad.pci,      pci];
            PSTH.events(n).bad.psth      = [PSTH.events(n).bad.psth,     betas];
            PSTH.events(n).bad.sem       = [PSTH.events(n).bad.sem,      sem];
            PSTH.events(n).bad.sessions  = [PSTH.events(n).bad.sessions, sessions(m)];
        end
    end
end
return
%------------------------------------------------------------------------

function [xY, SPM, SPMfile] = my_regions(xY,SPM)
% MY_REGIONS allows a user defined point or volume of interest. It returns
% a region structure containing the information about the region and its
% signal. NOTE: The entire time series is returned for later
% processing. It has been whitened and filtered.
%
% xY     - VOI structure (THIS IS SIMILAR BUT NOT IDENTICAL TO xY PRODUCED
%          BY SPM_REGIONS).
%       xY.xyz          - centre of VOI {mm}
%       xY.name         - name of VOI
%       xY.def          - VOI definition
%       xY.spec         - VOI definition parameters
%       xY.XYZmm        - MM coordinates of all voxels in VOI
%       xY.XYZ          - Voxel coordinates of all voxels in VOI
%       xY.y            - [whitened and filtered] voxel-wise data
%       xY.u            - first eigenvariate {scaled - c.f. mean response}
%                         SAME AS Y
%       xY.v            - first eigenimage
%       xY.s            - eigenvalues

if nargin < 1; xY = []; end;
if nargin < 2; SPM = []; end;
SPMfile = [];

if isempty(xY)

    xY.xyz  = spm_input('center of voi','!+1','e',[],3);
    xY.xyz  = xY.xyz(:);

    xY.name = spm_input('name of region','!+1','s','VOI');

    xY.def  = spm_input('VOI definition...','!+1','b',...
        {'point','sphere','box'});

    if strcmp(xY.def,'sphere') || strcmp(xY.def,'box')

        % assume we will need a result to get the voxels
        % unfortunately spm_getSPM is not scriptable
        % -------------------------------------------------------------
        spm_input('!SetNextPos',1);
        [SPM,xSPM] = spm_getSPM;
        Q          = ones(1,size(xSPM.XYZmm,2));
        SPMfile    = fullfile(SPM.swd,'SPM.mat');
    else
        if isempty(SPM)
            [SPM, SPMfile] = load_spm;
        end
    end

    % mm to voxel matrix
    % -------------------------------------------------------------
    iM = SPM.xVol.iM;

    switch xY.def
        case 'point'
            xY.spec = 0;
            xY.XYZ = iM * [xY.xyz; 1];
            xY.XYZ = xY.XYZ(1:3);
            xY.XYZmm = xY.xyz;

        case 'sphere'
            %---------------------------------------------------------------
            if ~isfield(xY,'spec')
                xY.spec = spm_input('VOI radius (mm)','!+0','r',0,1,[0,Inf]);
            end
            d        = [ xSPM.XYZmm(1,:) - xY.xyz(1);...
                xSPM.XYZmm(2,:) - xY.xyz(2);...
                xSPM.XYZmm(3,:) - xY.xyz(3) ];

            Q        = find(sum(d.^2) <= xY.spec^2);
            xY.XYZ   = xSPM.XYZ(:,Q);
            xY.XYZmm = xSPM.XYZmm(:,Q);

        case 'box'
            %---------------------------------------------------------------
            if ~isfield(xY,'spec')
                xY.spec = spm_input('box dimensions [x y z] {mm}',...
                    '!+0','r','0 0 0',3);
            end
            Q     = find(all(abs(xSPM.XYZmm - xY.xyz*Q) <= xY.spec(:)*Q/2));
            xY.XYZ = xSPM.XYZ(:,Q);
            xY.XYZmm = xSPM.XYZmm(:,Q);
    end
elseif isstruct(xY)

    % at this point we just assume that the user wants a point.
    % -------------------------------------------------------------
    if isfield(xY,'xyz')
        if isempty(SPM)
            [SPM, SPMfile] = load_spm;
        end

        xY.xyz  = xY.xyz(:);
        xY.spec = 0;
        xY.def  = 'point';
        if ~isfield(xY,'name')
            xY.name = 'VOI';
        end

        % mm to voxel matrix
        % -------------------------------------------------------------
        iM       = SPM.xVol.iM;
        xY.XYZ   = iM * [xY.xyz; 1];
        xY.XYZ   = xY.XYZ(1:3);
        xY.XYZmm = xY.xyz;
    end
else
    error('Not enough info to set up rest of xY structure');
end


%-Extract required data from results files
%=======================================================================

%-Get raw data, whiten and filter
%-----------------------------------------------------------------------
y = [];
try
    y   = spm_get_data(SPM.xY.VY,xY.XYZ);
    if isempty(y)

        % whoopsie no data
        % -------------------------------------------------------------
        error('No data in selected region.')
    end

    % whiten and filter
    % -------------------------------------------------------------
    y   = spm_filter(SPM.xX.K,SPM.xX.W*y);

catch
    try
        % remap files in SPM.xY.P if SPM.xY.VY is no longer valid
        %-------------------------------------------------------
        sprintf('Remapping data...')
        SPM.xY.VY = spm_vol(SPM.xY.P);
        y = spm_get_data(SPM.xY.VY,xY.XYZ);
        if isempty(y)

            % whoopsie no data
            % ----------------------------------------------------------
            error('No data in selected region.')
        end

        y = spm_filter(SPM.xX.K,SPM.xX.W*y);
    catch

        % data has been moved or renamed
        %-------------------------------------------------------
        y = [];
        spm('alert!',{'Original data have been moved or renamed',...
            'Recomendation: please update SPM.xY.P'},...
            mfilename,0);
    end
end

% compute regional response in terms of first eigenvariate
%-----------------------------------------------------------------------
[m n]   = size(y);
if m > n
    [v s v] = svd(spm_atranspa(y));
    s       = diag(s);
    v       = v(:,1);
    u       = y*v/sqrt(s(1));
else
    [u s u] = svd(spm_atranspa(y'));
    s       = diag(s);
    u       = u(:,1);
    v       = y'*u/sqrt(s(1));
end
d       = sign(sum(v));
u       = u*d;
v       = v*d;
Y       = u*sqrt(s(1)/n);

% set in structure
%-----------------------------------------------------------------------
xY.y    = y;
xY.u    = Y;
xY.v    = v;
xY.s    = s;

%========================================================================
function [SPM, SPMfile] = load_spm()

% Loads and SPM.mat fiile and returns the filename and the SPM structure
% itself.

SPM     = [];
SPMfile = [];

switch spm('ver')
    case 'SPM2'
        SPMfile = spm_get(1,'SPM.mat','Select SPM.mat');
    case 'SPM5'
        SPMfile = spm_select(1,'^SPM\.mat$','Select SPM.mat');
    otherwise
        error('Unknown version of SPM.');
end

[pth fn ex] = fileparts(SPMfile);
if ~isempty(pth)
    cd(pth)
else
    return
end
load(SPMfile);
spm_input('!SetNextPos',1);