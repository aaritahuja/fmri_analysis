% AA module - Iterative functional connectivity:
% John would like to start with MD candidate ROIs and use some sort of
% iterative optimization to improve upon them. Try to implement this, or
% something similar...
% djm 9/6/14
%
% 201114:
% - use Rik's updated script
%
% V2:
% - Now base threshold on RFX across all other animals (this may have broken some combinations of options?)
% - Now study level module, as will need to load all aubjects' data, but there are still per subject loops, so perhaps use spmd?
% - My rex.m updated to optionally combine svd threshold with max # components when choosing tissue confounds
% - roipath and roifilt can be paired cell arrays, in which case 1st cells give ROIs to be optimised,
%   and other cells give other sets of ROIs to include in s2s matrix in test mode
% - Don't apply bandpass and remove motion confounds before doing tissue confound svd, as it can make it unstable?
% - use single() for svd, because for unknown reason may fail to converge when double?
% - Remove options to do s2v only within GM mask, and to use canonical correlation
% TODO:
% - for test mode, don't need to load all subjects

function [aap,resp,zCssnames]=aamod_SPMconnIterateROIs2(aap,task,windowvec)
resp='';

switch task
    case 'domain'
        resp='study';
    case 'whentorun'
        resp='justonce';  % should this be run everytime or justonce?
    case 'description'
        resp='connectivity analysis';
    case 'summary'
        resp='connectivity analysis';
    case 'report'
        
    case {'doit','test'}
        
        % in test mode, just recreate SS matrix from left-out data of previously run iteration (need to recompute R as currently this is not saved)
        
        addpath /imaging/dm01/AndrewBell/caret_brain/djm/; % for GrowRegions3
        addpath /imaging/dm01/MoreTools/ % for fitplots2.m
        
        tissueconfound.history=aap.tasklist.currenttask.settings.tissueconfoundhistory;
        componentspertissue=aap.tasklist.currenttask.settings.componentspertissue;
        
        if strcmp(task,'doit')
            if ~isfield(aap.tasklist.currenttask.settings,'priors') ...
                    || isempty(aap.tasklist.currenttask.settings.priors)
                %priors='exponential'; % exponent would need to adapt to brain size; doesn't weight cluster size
                %priors='quadratic'; % falloff would need to adapt to brain size; doesn't weight cluster size
                %priors='smoothed'; % some function of ROI volume? (or brain size); slow to converge; poorly constrained?
                % may want to avoid or reduce overlap of priors?
                %priors='growthzone'; % hard constraint on range/direction of movement
                priors='adaptivegrowthzone';
                %priors='dilated'; % yet to try...
            else
                priors=aap.tasklist.currenttask.settings.priors;
            end
            if ismember(priors,{'exponential','quadratic','smoothed','growthzone','adaptivegrowthzone','dilated'})
                fprintf('\nUsing %s spatial priors',priors);
            elseif exist(priors,'file')
                fprintf('\nUsing "growth zones" constrained within "%s"',priors)
            else
                fprintf('\nUsing spatial priors from files suffixed "%s"',priors)
            end
        end
        
        if isempty(aap.tasklist.currenttask.settings.analname)
            analname='SPMconn';
        else
            analname=aap.tasklist.currenttask.settings.analname;
        end
        
        b=1;  % for now assume session 1
        
        %% get paths for all subjects
        subjpath=cell(1,length(aap.acq_details.subjects));
        sesspath=cell(1,length(aap.acq_details.subjects));
        outdir=cell(1,length(aap.acq_details.subjects));
        warn=true;
        for ws=1:length(aap.acq_details.subjects)
            subjpath{ws}=aas_getsubjpath(aap,ws);
            sesspath{ws}=aas_getsesspath(aap,ws,b);
            outdir{ws}=fullfile(subjpath{ws},analname);
            if strcmp(task,'doit')
                if ~exist(outdir{ws},'dir'),
                    mkdir(outdir{ws});
                else
                    if warn
                        uiwait(msgbox(sprintf('Contents of\n%s,\nand other subjects, will be deleted!',fullfile(outdir{ws})),'WARNING!','modal'));
                        warn=false;
                    end
                    delete(fullfile(outdir{ws},'*.*'));
                end
            end
        end
        groupdir=fullfile(aap.acq_details.root,analname);
        if ~exist(groupdir,'dir'),
            mkdir(groupdir);
        else
            delete(fullfile(groupdir,'*.*'));
        end
        
        %% prepare diary
        diaryname=fullfile(aap.acq_details.root,sprintf('diary_%s.txt',mfilename));
        if exist(diaryname,'file'), delete(diaryname); end
        diary(diaryname);
        
        %% get TR (may not be needed)
        try
            TR=aap.tasklist.currenttask.settings.TRs;
            if isempty(TR), gocatch; end
        catch
            try
                TR=aap.tasksettings.aamod_slicetiming.TRs; % TR (in seconds)
                if isempty(TR), gocatch; end
            catch
                warning('Assuming TR to be 2 s!')
                TR=2;
            end
        end
        
        %% get temporal filter and windows
        bpfilt=aap.tasklist.currenttask.settings.bpfilt;
        
        windows=aap.tasklist.currenttask.settings.windows;
        if windows==0 || ~isnumeric(windows), windows=1; end
        if ~exist('windowvec','var')
            windowvec=1:windows;
        end
        
        %% get epi prefix
        try
            if isfield(aap.tasklist.currenttask,'epiprefix') && ~isempty(aap.tasklist.currenttask.epiprefix)
                epiprefix=aap.tasklist.currenttask.epiprefix;
            else
                epiprefix=aap.tasklist.main.module(ismember({aap.tasklist.main.module.name},'aamod_SPMconnIterateROIs')).epiprefix;
            end
        catch
            warning('GUESSING EPIPREFIX TO BE "swar"')
            epiprefix='swar';
        end
        
        %% prepare matlabpool
        if matlabpool('size') ~= length(aap.acq_details.subjects),
            try
                matlabpool close;
            catch
                % sometimes it compains that one is open when there's not, 
                % perhaps if it crashes without closing properly??
            end
            fprintf('\n%s: ',datestr(now)); tic;
            try
                %%%%
                %             cwd=pwd;
                %             cd ~;
                %             myCluster = parcluster('CBU_Cluster');
                %             myCluster.NumWorkers=length(aap.acq_details.subjects);
                %             %myCluster.ResourceTemplate='-l nodes=^N^,mem=90GB,walltime=12:00:00';
                %             saveProfile(myCluster);
                %
                %             warning off all; matlabpool open; warning on all;
                %             cd(cwd);
                %%%%
                P=cbupool(length(aap.acq_details.subjects));
                P.ResourceTemplate='-l nodes=^N^,mem=360GB,walltime=12:00:00'; 
                % has completed 32 but not 35 subs with 90BG
                % is this total, or per lab?
                warning off all; matlabpool(P); warning on all;
                %%%%
                fprintf('Took %g mins',toc/60)
            catch err
                fprintf('Took %g mins, then failed to open matlabpool with error:',toc/60)
                error(err)
            end
        end
        
        txt={'GM','WM','CSF'};
        fprintf('\nLoading EPI data and masks')
        tic
        spmd % for ws=1:length(aap.acq_details.subjects)
            %% get epifiles and tissue masks for all subjects
            
            epifiles=aas_getimages(aap,labindex,b,epiprefix);
            if labindex==1, fprintf('.'); end
            allconfoundfiles=cellstr(epifiles);
            
            if ~isempty(aap.tasklist.currenttask.settings.confoundepiprefix)
                allconfoundfiles=regexprep(allconfoundfiles,...
                    sprintf('/%s([^/]*$)',epiprefix),...
                    sprintf('/%s$1',aap.tasklist.currenttask.settings.confoundepiprefix));
            end % might want to use different (e.g. smoothed) files for extracting confound (e.g. WM) timeseries
            if labindex==1, fprintf('.'); end
            epis=spm_vol(epifiles);
            if labindex==1, fprintf('.'); end
            
            switch lower(aap.tasklist.currenttask.settings.tissuepath)
                case 'uniqueindividual'
                    subid=regexprep(aap.acq_details.subjects(labindex).tag,'_.*','');
                    folder=fullfile(aap.acq_details.root,subid);
                otherwise
                    folder=fullfile(aas_getsubjpath(aap,labindex),'structurals');
            end
            GMmask=spm_select('FPlist',folder,aap.tasklist.currenttask.settings.GMfilt);
            WM=spm_select('FPlist',folder,aap.tasklist.currenttask.settings.WMfilt);
            CSF=spm_select('FPlist',folder,aap.tasklist.currenttask.settings.CSFfilt);
            %%%% for John, add extra ROIs from which to remove mean signal
            if ~isempty(aap.tasklist.currenttask.settings.extraconfoundroifilt)
                if exist([filesep aap.tasklist.currenttask.settings.extraconfoundroipath],'dir')
                    extrafolder=[filesep aap.tasklist.currenttask.settings.extraconfoundroipath];
                else
                    extrafolder=fullfile(aas_getsubjpath(aap,labindex),aap.tasklist.currenttask.settings.extraconfoundroipath);
                end
                extraconfoundroi=spm_select('FPlist',extrafolder,aap.tasklist.currenttask.settings.extraconfoundroifilt);
                warn=cellfun(@isempty,{GMmask,WM,CSF,extraconfoundroi});
            else
                warn=cellfun(@isempty,{GMmask,WM,CSF});
            end
            %%%%
            
            for w=find(warn)
                fprintf('\nWarning: subject %g: %s filter returned empty.',labindex,txt{w});
            end
            if labindex==1, fprintf('.'); end
             
            % analysis mask:
            if isfield(aap.tasklist.currenttask.settings,'analysismask') ...
                    && ~isempty(aap.tasklist.currenttask.settings.analysismask),
                mask=fullfile(outdir{labindex},'analysismask.nii');
                copyfile(aap.tasklist.currenttask.settings.analysismask,mask);
                spm_reslice(char(deblank(epifiles(1,:)),mask),struct('which',1,'mean',0,'interp',0,'prefix',''));
            else
                mask=spm_select('FPlist',fullfile(subjpath{labindex},'structurals'),'^wssbrainmask.nii$');
            end
        end
        fprintf(' (took %.2f mins)',toc/60)
        %         epifiles=epifiles(:);
        %         allconfoundfiles=allconfoundfiles(:);
        %         epis=epis(:);
        %         mask=mask(:);
        %         GMmask=GMmask(:);
        %         WM=WM(:);
        %         CSF=CSF(:);
        clear subid folder txt warn w
        
        %% seed ROIs:
        fprintf('\nLoading seed ROIs...'); tic
        spmd
            if isfield(aap.tasklist.currenttask.settings,'roipath') ...
                    && ~isempty(aap.tasklist.currenttask.settings.roipath)
                pth=aap.tasklist.currenttask.settings.roipath;
            else
                pth=fullfile(subjpath{labindex},'rois');
            end
            
            if ~iscell(pth),
                pth={pth};
                aap.tasklist.currenttask.settings.roifilt={aap.tasklist.currenttask.settings.roifilt};
            end % to allow multiple sets of ROIs for s2s matrix
            ROIs=cell(1,length(pth));
            for rs=1:length(pth)
                if ~exist([filesep,pth{rs}],'dir'), pth{rs}=fullfile(subjpath{labindex},pth); end
                ROIs{rs}=spm_select('FPlist',pth{rs},aap.tasklist.currenttask.settings.roifilt{rs}); %'^GM90.*.nii$'
            end
        end
        fprintf(' (took %.2f mins)',toc/60)
        %         aap=aap{1};
        %         ROIs=ROIs(:);
        
        %% don't need this until the iteration, but it's quite big so initialise it now
        if strcmp(task,'doit')
            spmd
                Yxyzrs=nan(epis(1).dim(1),epis(1).dim(2),epis(1).dim(3),size(ROIs{1},1),length(aap.acq_details.subjects),'single');
            end
        end
        
        %% get time windows and confounds
        %         allM=cell(1,length(aap.acq_details.subjects));
        %         scanind=cell(length(aap.acq_details.subjects),length(windowvec));
        %         R=cell(length(aap.acq_details.subjects),length(windowvec));
        
        fprintf('\nCreating confound regressors:')
        tic
        spmd
            
            % realignment parameters:
            allM=load(spm_select('FPlist',sesspath{labindex},'^rp.*.txt'));
            
            onsets=zeros(1,windows);
            durations=zeros(1,windows);
            scanind=cell(1,max(windowvec));
            for w=windowvec % windows
                %fprintf('\n\t=== Subject %g, Window %s',ws,num2str(w))
                if windows==1
                    onsets(1)=1;
                    durations(1)=length(allM)-1;
                else
                    onsets(w)=1+5+(w-1)*floor(length(allM)/windows); % note and consider gap between windows
                    durations(w)=floor(length(allM)/windows)-10;
                end
                scanind{w}=uint16(onsets(w):(onsets(w)+durations(w)));
                
                % prepare to call riks code to get adjusted motion/bandpass confounds,
                % in case we want to remove these before doing confound tissue svd
                S=[];
                S.StandardiseY=1;
                S.PreWhiten=0; % only if want stats
                S.pflag=0; % also partial regression?
                S.GlobalFlag=aap.tasklist.currenttask.settings.RemoveMeanGM;
                if aap.tasklist.currenttask.settings.globmove
                    % note allconfoundfiles is still on local client
                    S.GlobMove=y_FD_Jenkinson(allM(scanind{w},:),allconfoundfiles{1},40); % approx radius of brain = 40 mm
                    S.M=S.GlobMove;
                else S.M=M; %S.GlobMove;
                end
                S.VolterraLag=aap.tasklist.currenttask.settings.motionconfoundhistory;
                S.C=ones(length(scanind{w}),1); % just to get # scans for now...
                S.TR=TR;
                S.HPC=1/bpfilt(1);
                S.LPC=1/bpfilt(2);
                S.SpikeMovRelThr='';
                S.svd_thr=aap.tasklist.currenttask.settings.confoundsvdthresh;
                S.adjustonly=1;
                
                %[j, j, j, j, j, j, j, K, aM]=rsfMRI_GLM_djm(S,'quiet');
                
                %roi_threshold is hidden option!
                if isequal(componentspertissue,'mean')
                    [C1, j, params{1}, djm{1}]=rex(char(allconfoundfiles(scanind{w})),WM,'summary_measure','mean','roi_threshold',1-eps);
                    [C2, j, params{2}, djm{2}]=rex(char(allconfoundfiles(scanind{w})),CSF,'summary_measure','mean','roi_threshold',1-eps);
                else
                    % cant find any explanation of whether to prefer pca=0 or =1. 1 is default, but conn seems to use 0
                    % "PCA decomposition 1 = first extract component represents mean signal, next components represent principal component decomposition of covariance matrix);
                    % 0 = components represent SVD decomposition of the second moment (X'*X) matrix"
                    % occasionally fails with covariates [K aM] (and 'pca',0), but not [aM K] ?!
                    [C1, j, params{1}, djm{1}]=rex(char(allconfoundfiles(scanind{w})),WM,'summary_measure','eigenvariate','dims',componentspertissue,'pca',0,'roi_threshold',1-eps,'covariates',[]);
                    [C2, j, params{2}, djm{2}]=rex(char(allconfoundfiles(scanind{w})),CSF,'summary_measure','eigenvariate','dims',componentspertissue,'pca',0,'roi_threshold',1-eps,'covariates',[]);
                end
                
                %%% djm: save tissue confound basis images
                nrow=min(10,max(size(params{1}.ROIdata,2), size(params{2}.ROIdata,2)));
                fig=spm_figure('getwin'); clf
                set(fig,'units','normalized');
                pos=fitplots2([nrow 2]);
                for tiss=1:2
                    XYZ=params{tiss}.ROIinfo.voxels{1}{1};
                    for c=1:min(nrow,size(params{tiss}.ROIdata,2))
                        Z=params{tiss}.ROIinfo.basis{1}{1}(:,c);
                        mip=load('/imaging/dm01/AndrewBell/RhesusMacaqueAtlasTemplates/atlasfiles/RhesusMIP.mat');
                        if isfield(mip,'mask_all'), mip=1+mip.mask_all; else mip=1+mip.mip96; end
                        if length(unique(Z))==1, Z=2+zeros(size(Z));
                        else Z=2+Z/max(abs(Z)); end
                        d=spm_project(Z',round(XYZ'),[2,2,2,size(mip)]);
                        idx=find(d~=0);mip(idx)=round(34.5+31.5*(d(idx)-2));
                        mip(all(mip==2,2),:)=[];
                        mip(:,all(mip==2,1))=[];
                        subplot('position',pos{c,tiss})
                        image(rot90(mip));axis equal;axis tight;axis off;
                    end % next eigenvector
                end % next tissue class ROI
                colormap([1*ones(1,3);.8*ones(1,3);jet(64)]);
                if isequal(componentspertissue,'mean')
                    txt='MIPs of WM and CSF components';
                else
                    txt=sprintf('MIPs of top %g (of %g, expaining %.1f%% var) WM,  and top %g (of %g, expaining %.1f%% var) CSF components',...
                        min(nrow,size(params{1}.ROIdata,2)),size(params{1}.ROIdata,2),...
                        djm{1}.varexplained*100,...
                        min(nrow,size(params{2}.ROIdata,2)),size(params{2}.ROIdata,2),...
                        djm{2}.varexplained*100);
                end
                tmp=annotation('textbox',[0 0.95 1 0.05],'Color','k','String',txt,'LineStyle','none');
                print(fig, fullfile(outdir{labindex},sprintf('WMCSFcomponents_window%g.png',w)),'-dpng');
                delete(tmp)
          
                if tissueconfound.history==1 % first-order difference
                    aC1=[C1, [zeros(1,size(C1,2)); diff(C1)]];
                    aC2=[C2, [zeros(1,size(C2,2)); diff(C2)]];
                elseif tissueconfound.history>1 % Volterra expansion as Rik uses for motion (but with initial padding)
                    % gives 90 columns, but rank ~20 - correct?
                    bf = eye(tissueconfound.history);  % artifacts can last up to 5 TRs = 10s, according to Power et al (2013)
                    bf = [bf; diff(bf)];
                    U=[]; for c=1:size(C1,2); U(c).u = C1([ones(1,tissueconfound.history), 1:size(C1,1)],c); U(c).name{1}='c'; end;
                    aC1 = spm_Volterra(U,bf',2);
                    aC1(1:tissueconfound.history,:)=[];
                    U=[]; for c=1:size(C2,2); U(c).u = C2([ones(1,tissueconfound.history), 1:size(C2,1)],c); U(c).name{1}='c'; end;
                    aC2 = spm_Volterra(U,bf',2);
                    aC2(1:tissueconfound.history,:)=[];
                end
                %aC = spm_en([aC1 aC2],0);
                aC =[aC1 aC2];
                
                fprintf('\n\t expanded WM/CSF confounds using %s dfs; ',mat2str([size(aC1,2) size(aC2,2)]));
                
                %%%% for John, add extra ROIs from which to remove mean signal
                if ~isempty(aap.tasklist.currenttask.settings.extraconfoundroifilt)
                    fprintf('\n\t adding mean signal from %s as extra confound',extraconfoundroi)
                    aC=[aC rex(char(allconfoundfiles(scanind{w})),extraconfoundroi,'summary_measure','mean','roi_threshold',1-eps)];
                end
                
                % call riks code again to get adjusted confounds after including tissue confounds
                S.C=aC;
                [j, j, j, j, j, X0s]=rsfMRI_GLM_djm(S);
                X0s=X0s(:,1:end-1); % remove constant? Identical correlations
                
                try
                    R{w}  = eye(length(scanind{w})) - X0s*pinv(single(X0s));
                catch
                    warning('Using svds for pinv')
                    R{w}  = eye(length(scanind{w})) - X0s*pinv_svds(X0s);
                end
                
            end % next window
        end
        fprintf('Creating confound regressors took %.2f mins.',toc/60)
                 clear onsets durations S *C1 *C2 aC j params djm nrow fig pos XYZ c Z mip d idx txt tmp bf U X0s w
                 clear allM tiss
        %         R=R(:);
        %         scanind=scanind(:);
        
        %% extract/summarise data from all voxels of all subjects
        % (this is here so we can load data once then process multiple time windows, but might hit memory limits?)
        if strcmp(task,'doit')
            fprintf('\nExtracting all brain voxel and gm timeseries...');
            tic;
            spmd
                %if ~strcmp(task,'doit') && ws~=i, continue; end
                [allvox, allXYZvox]=getdata(epis,mask,'','','',labindex>1 || strcmp(task,'test')); % whole brain mask
                %[gmvox, gmXYZvox]=getdata(epis,GMmask,'','','',labindex>1 || strcmp(task,'test'));
                %         allvox=allvox(:);
                %         allXYZvox=allXYZvox(:);
                %         gmvox=gmvox(:);
                %         gmXYZvox=gmXYZvox(:);
            end
            fprintf(' (took %.2f mins)',toc/60)
        end
        
        %% economise ROI names
        % for main iterated seeds:
        tmp=ROIs{1}; % Composite objects only support simple subscripting.
        nroisets=length(tmp);
        ROIlabels=cell(1,nroisets);
        bits=cell(1,nroisets);
        for rs=1:nroisets
            [junk, bits{rs}]=spm_str_manip(tmp{rs},'C');
            ROIlabels{rs}=bits{rs}.m;
        end
        % L/R and extra ROIs for s2s matrix:
        zCssnames={};
        for rs=1:length(ROIlabels)
            zCssnames=[zCssnames; strcat('L_',ROIlabels{rs}); strcat('R_',ROIlabels{rs})];
        end
        
        %% copy ROIs to output folder (will assume they are binary) and reslice to EPI
        % (do for all ROI sets, even though only 1st ROI set will be optimised,
        % because want them treated equally in s2s matrix i.e. all resliced in same way)
        fprintf('\nCopying/reslicing initial ROIs...');
        tic
        spmd
            
            roifilesO=cell(1,nroisets);
            for rs=1:nroisets
                roifilesO{rs}=cell(1,length(ROIlabels{rs}));
                for r=1:length(ROIlabels{rs})
                    if strcmp(task,'doit')
                        roifilesO{rs}{r}=fullfile(outdir{labindex},sprintf('%s_iteration0.nii',ROIlabels{rs}{r}));
                    else
                        roifilesO{rs}{r}=fullfile(outdir{labindex},sprintf('%s_TEST%s.nii',ROIlabels{rs}{r},aap.tasklist.currenttask.settings.doneflagsuffix));
                    end
                    copyfile([bits{rs}.s ROIlabels{rs}{r} bits{rs}.e],roifilesO{rs}{r});
                    if labindex==1, fprintf('.'); end
                end
                %ROIlabels{rs}=strcat('ROI_',ROIlabels{rs}); % prefix in case bits start with a number and will become field names
                spm_reslice(char(deblank(epifiles(1,:)),char(roifilesO{rs})),struct('which',1,'mean',0,'interp',0,'prefix',''));
            end
        end % next subject
        fprintf(' (took %.2f mins)',toc/60)
        
        %% combine initial ROIs
        if strcmp(task,'doit')
            fprintf('\nCombining initial ROIs...');
            tic
            spmd
                %%% put main set into single nii to compare across iterations
                Vout=rmfield(epis(1),'pinfo');
                MPPO=zeros(Vout.dim);
                vi=sub2ind(Vout.dim,allXYZvox(1,any(allvox)),allXYZvox(2,any(allvox)),allXYZvox(3,any(allvox)));
                MPPO(vi)=1;
                
                for r=1:length(ROIlabels{1})
                    temp=spm_read_vols(spm_vol(roifilesO{1}{r}));
                    MPPO(vi)= MPPO(vi) + r*(temp(vi) > 0);
                    if labindex==1, fprintf('.'); end
                end
                if strcmp(task,'doit')
                    MPPfile=fullfile(outdir{labindex},'MPP_iteration0.nii');
                else
                    MPPfile=fullfile(outdir{labindex},sprintf('MPP_TEST%s.nii',aap.tasklist.currenttask.settings.doneflagsuffix));
                end
                Vout.fname=MPPfile;
                Vout.dt=[4 0];
                spm_write_vol(Vout,MPPO);
                MPPO=uint8(MPPO);
                %labBarrier; % ensure that all files get written before they are accessed
            end
            fprintf(' (took %.2f mins)',toc/60)
        end
        
        %% define the spatial priors (currently assuming constant across subjects)
        if strcmp(task,'doit')
            fprintf('\nDefining spatial priors (assuming constant across subjects):');
            tmppth=pth{1};  % Composite objects only support simple subscripting.
            tmpVout=Vout{1};
            aap=aap{1};
            if exist(priors,'file')
                incmask=priors;
                priors='growthzone';
            elseif ~ismember(priors,{'exponential','quadratic','smoothed','growthzone','adaptivegrowthzone','dilated'})
                priors=cellstr(spm_select('FPlist',tmppth{1},regexprep(aap.tasklist.currenttask.settings.roifilt{1},'(.nii|.img)',[priors '$1'])));
            end
            
            halo=cell(1,length(ROIlabels{1}));
            if iscell(priors)
                catchmentO=zeros(tmpVout.dim);
                overlap=[];
                weight=zeros(length(priors),numel(catchmentO));
                for p=1:length(priors)
                    [j, co,j, co0]=getdata(MPPfile{1},priors{p}); % find coordinates of voxels where prior is in brain mask
                    hit=sub2ind(tmpVout.dim,co(1,:),co(2,:),co(3,:));
                    for c=1:length(hit)
                        weight(p,hit(c))=sum(ismember(round(co0'),co(:,c)','rows'));
                    end
                    check=hit(catchmentO(hit)>0);
                    overlap=[overlap, check]; % to fix voxels that overlap >1 prior
                    catchmentO(hit)=p+1;
                end
                [j,ind]=max(weight);
                catchmentO(overlap)=ind(overlap)+1;
                tmpVout.fname=fullfile(groupdir,'LabelledPriorRegions.nii');
                spm_write_vol(tmpVout,catchmentO);
                % write out individual growth zones (might be useful later)
                maskname=cell(1,length(ROIlabels{1}));
                for r=1:length(ROIlabels{1})
                    maskname{r}=regexprep(fullfile(groupdir,'LabelledPriorRegions.nii'),'....$',sprintf('_%s_only.nii',ROIlabels{1}{r}));
                    tmpVout.fname=maskname{r};
                    spm_write_vol(tmpVout,catchmentO==(r+1));
                end
                % write out 'halo' of individual growth zones, in case want to test connectivity within ROI against a baseline that's more local than e.g. whole brain GM
                for r=1:length(ROIlabels{1})
                    halofn=regexprep(fullfile(groupdir,'LabelledPriorRegions.nii'),'....$',sprintf('_%s_only_halo.nii',ROIlabels{1}{r}));
                    tmpVout.fname=halofn;
                    halo{r}=logical(spm_dilate(double(catchmentO==(r+1)))-(catchmentO==(r+1)));
                    %spm_write_vol(tmpVout,halo{r});
                end
                priors='masks';
            else
                switch lower(priors)
                    case {'growthzone','adaptivegrowthzone'}
                        % growth of MPP (potential fixed prior)
                        clear temp
                        temps.labelimage=true;
                        temps.coresubregions{1}{1}=MPPfile{1};
                        temps.ignore=1; % analysis mask
                        temps.minvox=1;
                        temps.anticore{1}='';
                        if exist('incmask','var') && exist(incmask,'file')
                            temps.maskI{1}=incmask;
                        else
                            temps.maskI{1}=mask; % whole brain mask
                        end
                        temps.maskE{1}='';
                        temps.maxgrowth=Inf;
                        temps.maxclusters=2;
                        temps.pref='Grown';
                        V=GrowRegions3(temps);
                        catchmentO=round(spm_read_vols(V));
                        % write out individual growth zones (might be useful later)
                        tmpVout=rmfield(V,'pinfo');
                        for r=1:length(ROIlabels{1})
                            tmpVout.fname=regexprep(V.fname,'....$',sprintf('_%s_only.nii',ROIlabels{1}{r}));
                            spm_write_vol(tmpVout,catchmentO==(r+1));
                        end
                        % write out 'halo' of individual growth zones, in case want to test connectivity within ROI against a baseline that's more local than e.g. whole brain GM
                        for r=1:length(ROIlabels{1})
                            halofn=regexprep(V.fname,'....$',sprintf('_%s_only_halo.nii',ROIlabels{1}{r}));
                            tmpVout.fname=halofn;
                            halo{r}=logical(spm_dilate(double(catchmentO==(r+1)))-(catchmentO==(r+1)));
                            %spm_write_vol(tmpVout,halo{r});
                        end
                end % which prior type?
            end % priors is cell array?
            clear c r V P j co co0 check hit overlap incmask weight ind p
            clear CSF WM X0s allconfoundfiles epifiles halofn junk pth rs sesspath subjpath temp tmpVout tmppth vi
            %% make size counter for each ROI in each subject
            spmd
                for n=1:length(ROIlabels{1})
                    nvox(n)=sum(round(MPPO(:))==(n+1));
                end
            end
        end % do it?
        
        %% Adjust whole-brain data for confounds
        if strcmp(task,'doit')
        fprintf('\nAdjusting data for confounds...');
        tic
        spmd
            for w=windowvec(:)'
                %%% create adjusted data, all voxels, all subjects:
                aYvox{w}=R{w}*allvox(scanind{w},:);
            end
        end
        fprintf(' (took %.2f mins)',toc/60)
        end
        
        %% collect or codistribute variables so they are available to all workers
        fprintf('\nCombining/redistributing variables from workers...') % because all workers/subjects now need access to data from all other workers/subjects
        tic
        %%%
        %         epis={epis{:}}; % clumsy
        %         scanind={scanind{:}};
        %         roifilesO={roifilesO{:}};
        %         R={R{:}};
        %         allvox={allvox{:}};
        %         allXYZvox={allXYZvox{:}};
        %%%
        %         spmd
        %             xepis=codistributed.build({epis},codistributor('1d',2,ones(1,length(aap.acq_details.subjects)),[1,length(aap.acq_details.subjects)]));
        %             xscanind=codistributed.build({scanind},codistributor('1d',2,ones(1,length(aap.acq_details.subjects)),[1,length(aap.acq_details.subjects)]));
        %             xroifilesO=codistributed.build({roifilesO},codistributor('1d',2,ones(1,length(aap.acq_details.subjects)),[1,length(aap.acq_details.subjects)]));
        %             xroifiles=codistributed.build({roifilesO},codistributor('1d',2,ones(1,length(aap.acq_details.subjects)),[1,length(aap.acq_details.subjects)]));
        %             xR=codistributed.build({R},codistributor('1d',2,ones(1,length(aap.acq_details.subjects)),[1,length(aap.acq_details.subjects)]));
        %             xallvox=codistributed.build({allvox},codistributor('1d',2,ones(1,length(aap.acq_details.subjects)),[1,length(aap.acq_details.subjects)]));
        %             xallXYZvox=codistributed.build({allXYZvox},codistributor('1d',2,ones(1,length(aap.acq_details.subjects)),[1,length(aap.acq_details.subjects)]));
        %             xgmvox=codistributed.build({gmvox},codistributor('1d',2,ones(1,length(aap.acq_details.subjects)),[1,length(aap.acq_details.subjects)]));
        %             xgmXYZvox=codistributed.build({gmXYZvox},codistributor('1d',2,ones(1,length(aap.acq_details.subjects)),[1,length(aap.acq_details.subjects)]));
        %             xGMmask=codistributed.build({GMmask},codistributor('1d',2,ones(1,length(aap.acq_details.subjects)),[1,length(aap.acq_details.subjects)]));
        %         end
        %         clear epis scanind roifilesO R allvox allXYZvox GMmask gmvox gmXYZvox
        %%%
        roifilesO=roifilesO(:);
        R=R(:);
        GMmask=GMmask(:);
        
        if strcmp(task,'doit')
            epis=epis(:);
            % for n==35, need to save space:
            epis=cellfun(@(x) {x.fname},epis,'uniformoutput',0);
            scanind=scanind(:);
            allXYZvox=allXYZvox(:);
            %gmvox=gmvox(:);
            %gmXYZvox=gmXYZvox(:);
            
            aYvox=aYvox(:);
            spmd
                aYvox=aYvox;
            end
        end
        %%%
        fprintf(' (took %.2f mins)',toc/60)
        clear gmvox gmXYZvox allvox*
        
        %% main iteration per window:
        for w=windowvec(:)' %windows
            
            epinames=cell(size(epis));
            for e=1:length(epis) % for n==35, need to save space
                epinames{e}=epis{e}(scanind{e}{w});
            end
            
            if strcmp(task,'doit')
                fprintf('\n == Time window %g ==',w);
                switch lower(priors)
                    case {'growthzone','adaptivegrowthzone','masks'}
                        catchment=catchmentO;
                end
            end
            
            fprintf('\nRunning iterations for % g subjects and %g seeds. Started %s\n', length(aap.acq_details.subjects),length(ROIlabels{1}),datestr(now))
           tic
%%%
%             MPPO=MPPO{1};
%             MPPfile=MPPfile{1};
%             ROI=ROIs{1};
%             Vout=Vout{1};
%             mask=mask{1};
%             n=n{1};
%             nvox=nvox{1};
%             matlabpool close
%             %%%
            spmd % for ws=1:length(aap.acq_details.subjects)
                % for each subject, run iteration using all other subjects:
                
%                 % prepare outputs (to be used per iteration)
%                 outvol(1:length(ROIlabels{1}),1:length(aap.acq_details.subjects))=rmfield(epis{labindex}(1),'pinfo');
%                 [outvol.fname]=deal('');
%                 [outvol.dt]=deal([16 0]); % to keep nans
                %outvolnames=cell(1:length(ROIlabels{1}),1:length(aap.acq_details.subjects));
                
                iterate=true; it=0; ham=0;
                tobreak=0;
                while iterate
%                     if all(gcat(tobreak)),
%                         % for version A. Does this work?
%                         break; % get out of iteration loop only when all labs have finished
%                     end
                    
                    if tobreak,
                        % this lab has finished, but others haven't.
                        % For version A, don't do anything, but stay in loop so labbarrier can work.
                        % For version B, now is the time to break
                        break
                    else
                        it=it+1;
                        %if labindex==1
                            fprintf('\nIteration %03d (taken %g mins). Getting connectivity for all subs',it,round(toc/60))
                        %end
                        %%% get connectivity to all voxels from each source ROI, from all other subjects
                        %Yxyzrs=nan(epis{labindex}(1).dim(1),epis{labindex}(1).dim(2),epis{labindex}(1).dim(3),length(ROIlabels{1}),length(aap.acq_details.subjects),'single');
                        Yxyzrs(:)=nan;
                        
                        %%%%%% Main ROIs Version B: Iterated ROIs from current subject are always applied to all other subjects,
                        % to keep current subject isolated across all iteration steps
                        if it==1
                            lastfiles=roifilesO{labindex}{1};
                        else
                            lastfiles=regexprep(roifilesO{labindex}{1},'_iteration0',sprintf('_window%g_iteration%g',w,it-1));
                        end
                        %%%%%%
                        
                        for wis=setdiff(1:length(aap.acq_details.subjects),labindex) % each other subject
                            %%% create adjusted data
                            
                            %                         %%%%% main ROIs: Version A: Iterated ROIs from each other subject are used.
                            %                         % This means that although data from current subject are held out on each iteration,
                            %                         % the data are mixed over multiple iteration steps. Also it might currently fail because
                            %                         % it tries to load later iteration steps from a subject whose iteration has terminated
                            %                         if it==1
                            %                             lastfiles=roifilesO{wis}{1};
                            %                         else
                            %                             lastfiles=regexprep(roifilesO{wis}{1},'_iteration0',sprintf('_window%g_iteration%g',w,it-1));
                            %                         end
                            %                         %%%%%%
                            
                            aY=getdata(epinames{wis},lastfiles,'nanmean',[],GMmask{wis},true);
                            if labindex==1,fprintf('.'); end
                            aY = R{wis}{w}*aY; %figure(16); imagesc([zscore(Y)'; aY']); colorbar
                            
                            %%% seed to voxel
                            for n=1:length(ROIlabels{1}) % each seed
                                zC=corr(aY(:,n),aYvox{wis}{w});
                                zC=atanh(zC); % fisher transform
                                
                                %outvol(n,wis).fname=fullfile(outdir{labindex},sprintf('TEMP_%s_window%g_zC_sub%g.nii',ROIlabels{1}{n},w,wis));
                                %outvolnames{n,wis}=fullfile(outdir{labindex},sprintf('TEMP_%s_window%g_zC_sub%g.nii',ROIlabels{1}{n},w,wis));
                                
                                ind=sub2ind(Vout.dim,allXYZvox{wis}(1,:),allXYZvox{wis}(2,:),allXYZvox{wis}(3,:));
                                Yout=nan(Vout.dim);
                                Yout(ind)=zC;
                                Yxyzrs(:,:,:,n,wis)=single(Yout);
                            end % next main seed
                            if labindex==1, fprintf('\b:'); end
                            aY=[]; zC=[]; ind=[]; 
                        end % next subject
                        
                        %%% redefine each ROI
                        if labindex==1, fprintf(' Redefining ROIs'); end
                        Nr=length(ROIlabels{1});
                        PP=nan([size(Yout) Nr],'single'); % the error: "The variable contains a hidden composite..." does not arise due to how the variable is initialised, but how it might be used...
                        Yout=[];
                        for n=1:Nr
                            %V=outvol(n);
                            
                            %%% average voxelwise connectivity maps from other ROIs (and rescale?)
                            Y=nanmean(Yxyzrs(:,:,:,setdiff(1:Nr,n),:),4);
                            Y(Y==0)=nan; % outside brain mask
                            Y=squeeze(tanh(Y));
                            %%% get mean connectivity across non-ROI voxels from other ROIs and other subjects
                            temp=nan(1,length(aap.acq_details.subjects));
                            for wis=setdiff(1:length(aap.acq_details.subjects),labindex) % each other subject
                                %%% all gm?:
                                tempY=Y(:,:,:,wis);
                                %ind=sub2ind(outvol(n,wis).dim,gmXYZvox{wis}(1,:),gmXYZvox{wis}(2,:),gmXYZvox{wis}(3,:));
                                %temp(wis)=nanmean(tempY(ind));
                                %%% or halo of current ROI?:
                                temp(wis)=nanmean(tempY(halo{n}>0));
                                tempY=[];
                            end
                            baseline=nanmean(temp);
                            %%% do RFX across subjects to test for higher connectivity than this mean connectivity value
                            [h, p]=ttest(Y,baseline,'alpha',0.05,'dim',4,'tail','right');
                            Y=fdr_bh(p,0.05,'pdep');
                            %%%
                            
                            switch priors
                                case {'exponential','quadratic'}
                                    Y=Y/nansum(Y(:))*nvox; % and scale to sum to nvox within brainmask
                            end
                            
                            %%% make spatial prior (from distance transform, or smoothing,
                            % or growth zone of candidate)
                            switch priors
                                case 'exponential'
                                    if ~nvox(it,n) % use last available prior:
                                        D=spm_read_vols(spm_vol(fullfile(outdir,sprintf('%s_window%g_iteration%g_ROIdistprior.nii',ROIlabels{1}{n},w,it-1))));
                                    else
                                    D=1./(bwdist(R{n})); % exponential
                                    D(isinf(D))=1;
                                    D(isnan(Y))=nan;
                                    D=D/nansum(D(:))*nvox; % and scale to sum to nvox within brainmask
                                    end
                                case 'quadratic'
                                    if ~nvox(it,n) % use last available prior:
                                        D=spm_read_vols(spm_vol(fullfile(outdir,sprintf('%s_window%g_iteration%g_ROIdistprior.nii',ROIlabels{1}{n},w,it-1))));
                                    else
                                    D=bwdist(R{n}); % linear...
                                    D(isnan(Y))=nan;
                                    D=nanmax(D(:))-D; % flip to closeness
                                    D=D.^2; % ...or quadratic
                                    D=D/nansum(D(:))*nvox; % and scale to sum to nvox within brainmask
                                    end
                                case 'smoothed'
                                    if ~nvox(it,n) % use last available prior:
                                        D=spm_read_vols(spm_vol(fullfile(outdir,sprintf('%s_window%g_iteration%g_ROIdistprior.nii',ROIlabels{1}{n},w,it-1))));
                                    else
                                    D=zeros(size(R{n})); % smoothing
                                    spm_smooth(double(R{n}),D,diam/4)
                                    D=D/nanmax(D(:)); % and scale so max is 1
                                    end
                                case {'growthzone','adaptivegrowthzone','masks'}
                                    D=(catchment==(n+1)); % +1 for analysis mask
                                case 'dilated'
                                    %%% not implemented
                            end
                            
                            %%% combine
                            switch priors
                                case {'exponential','quadratic'}
                                    PP(:,:,:,n)=Y/nvox(n,it) .* D/nvox(n,it);
                                otherwise
                                    PP(:,:,:,n)=Y .* D;
                            end
                            if labindex==1, fprintf('.'); end
                            D=[];
                        end % next ROI
                        
                        switch lower(priors)
                            case 'adaptivegrowthzone'
                                % regrowth of MPP
                                opt=struct([]);
                                opt.labelimage=true;
                                opt.coresubregions{1}{1}=MPPfile;
                                opt.ignore=1; % analysis mask
                                opt.minvox=1;
                                opt.anticore{1}='';
                                opt.maskI{1}=mask; % whole brain mask
                                opt.maskE{1}='';
                                opt.maxgrowth=Inf;
                                opt.maxclusters=Inf;
                                opt.pref='Grown';
                                tempV=GrowRegions3(opt);
                                catchment=round(spm_read_vols(tempV));
                        end
                        
                        %if it==1, O=R; end % store (O)riginal (R)OIs, to measure subsequent change
                        
                        [maxP, MPP]=max(PP,[],4); % maxP will probably be positive because multiplying with the prior will tend to put things to zero
                        PP=[];
                        
                        %%% create the combined ROI file
                        MPPfile=fullfile(outdir{labindex},sprintf('MPP_window%g_iteration%g.nii',w,it));
                        V=Vout;
                        V.fname=MPPfile;
                        V.dt=[4 0];
                        MPP=MPP .* ~isnan(maxP) .* ~(maxP==0) + ~isnan(Y);
                        % zero outside brain, then zero outside ROIs, then +1 in brain or GM so "non-network" will be extracted by getdata
                        spm_write_vol(V, MPP );
                        maxP=[]; Y=[];
                        
                        % save individual volumes & plot ROI size
                        figure(11); clf
                        positions=fitplots2(Nr+1);
                        for n=1:Nr
                            roifilesO{labindex}{1}{n}=fullfile(outdir{labindex},sprintf('%s_window%g_iteration%g.nii',ROIlabels{1}{n},w,it));
                            V.fname=roifilesO{labindex}{1}{n};
                            spm_write_vol(V,MPP==(n+1));
                            nvox(it+1,n)=nansum(MPP(:)==(n+1));
                            
                            subplot('position',positions{n});
                            sizechange=round( (nvox(:,n)-nvox(1,n)) /nvox(1,n)*100 );
                            jd(it,n)=pdist(double([MPP(:)==(n+1),MPPO(:)==(n+1)]'),'jaccard');
                            if it>1
                                ax=plotyy(0:it,sizechange',0:it,[0 jd(:,n)']);
                                ylabel(ROIlabels{1}{n},'interpreter','none')
                                set(ax,'xtick',[]);
                                %axis tight; ylim([-100 200])
                            end
                        end % next roi
                        if it>1
                            set(ax,'xtick',0:it,'xticklabel',0:it)
                            legend({'ROI size (% change)','Jaccard distance'},'location','nw')
                        end
                        
                        % calculate change from (any) previous set
                        tobreak=false;
                        MPP=spm_read_vols(spm_vol(fullfile(outdir{labindex},sprintf('MPP_window%g_iteration%g.nii',w,it))));
                        % not sure why I have to reload this, but something seems to change when saving and loading
                        % (perhaps at edge of volume??) which would otherwise prevent iteration from terminating
                        for previous=(it-1):-1:0
                            if ~previous
                                %b=MPPO;
                                b=spm_read_vols(spm_vol(fullfile(outdir{labindex},sprintf('MPP_iteration0.nii'))));
                            else
                                b=spm_read_vols(spm_vol(fullfile(outdir{labindex},sprintf('MPP_window%g_iteration%g.nii',w,previous))));
                            end
                            good=(MPP(:)>0 & b(:)>0);
                            check=pdist([MPP(good), b(good)]','hamming'); % identical to pdist([a(:), b(:)]','jaccard')
                            % described as percentage but actually seems to be proportion
                            if ~check, tobreak=true; end
                        end
                        % have got stuck in a loop where e.g. two voxels oscillate
                        
                        ham(it)=check;
                        if it>1
                            subplot('position',positions{end});
                            plot([0 ham],'-k.','color',[0 .6 0]); axis tight
                            ylabel('Overall')
                            set(gca,'xtick',0:it,'xticklabel',0:it);
                            legend({'Hamming distance'},'location','se')
                        end
                        
                        if labindex==1
                            fprintf(' sizes: %s  Hamming dist from original=%.3f',mat2str(nvox(it,:)),check)
                            %                         if it>1
                            %                             a=spm_read_vols(spm_vol(fullfile(outdir{labindex},sprintf('MPP_window%g_iteration%g.nii',w,it))));
                            %                             b=spm_read_vols(spm_vol(fullfile(outdir{labindex},sprintf('MPP_window%g_iteration%g.nii',w,it-1))));
                            %                             check2=pdist([a(good), b(good)]','hamming');
                            %                             fprintf('\nHamming distance from previous after reloading=%.3f',check2)
                            %                         end
                        end
                        
                    end % do stuff if this lab hasn't stopped iteration
                    
                    % needed for version A:?
                    %labBarrier % force all labs to catch up before next iteration
                    
                end % iterating
                
                print(11,'-dpng',fullfile(outdir{labindex},sprintf('Iterations_window%g.png',w)));
                
            end % next subject
            fprintf(' (took %.2f mins)',toc/60)
            
            %%%%%%%%%% save seed to voxel maps from first and last iteration of test subjects
            fprintf('\nSaving seed-to-voxel maps for original and optimised ROIs, and masks:');
            tic;
            clear lastfiles fn
            spmd
                %lastfiles=cell(2,length(roifilesO{labindex}{1}));
                %fn=cell(1,length(ROIlabels));
                for wit=[1 2 3]; % first/mask/last
                    if wit==1
                        lastfiles{1,:}=regexprep(roifilesO{labindex}{1},'_window\d+_iteration\d+.nii','_iteration0.nii');
                    elseif wit==2
                        lastfiles{2,:}=maskname;
                    else
                        lastfiles{3,:}=roifilesO{labindex}{1};
                    end
                    
                    temproidata=getdata(epinames{labindex},lastfiles{end,:},'nanmean',[],GMmask{labindex},true);
                    aY = R{labindex}{w}*temproidata; %figure(16); imagesc([zscore(Y)'; aY']); colorbar
                    
                    %%% seed to voxel
                    for rs=1 %:length(ROIlabels) % just from seeds, not from control regions
                        for n=1:length(ROIlabels{rs}) % each seed
                            C=corr(aY(:,n),aYvox{labindex}{w});
                            zC=atanh(C); % fisher transform
                            if wit==2
                                fn{rs}{2,n}=fullfile(outdir{labindex},sprintf('%s_window%g_zC_mask.nii',ROIlabels{rs}{n},w));
                            elseif wit==1
                                fn{rs}{wit,n}=fullfile(outdir{labindex},sprintf('%s_window%g_zC_iteration0.nii',ROIlabels{rs}{n},w));
                            else
                                fn{rs}{wit,n}=fullfile(outdir{labindex},sprintf('%s_window%g_zC_iteration%g.nii',ROIlabels{rs}{n},w,it));
                            end
                            V.fname=fn{rs}{wit,n};
                            
                            ind=sub2ind(Vout.dim,allXYZvox{labindex}(1,:),allXYZvox{labindex}(2,:),allXYZvox{labindex}(3,:));
                            Yout=nan(Vout.dim);
                            Yout(ind)=zC;
                            spm_write_vol(V,Yout);
                            if labindex==1, fprintf('\b:'); end
                        end % next seed
                    end % nest seed set
                end % first and last iterations
            end % next subject
            fprintf(' (took %.2f mins)',toc/60)
            
            %% save filenames for ease of group analysis outside of this function
            connectivityfiles=fn(:);
            ROIfiles=lastfiles(:);
            save(fullfile(groupdir,'FinalFiles.mat'),'connectivityfiles','ROIfiles')
            
            %% create mean connectivity maps across subjects, and do RFX
            fprintf('\nCreating mean connectivity maps')
            RMSL_to_F99_job='/imaging/dm01/AndrewBell/caret_brain/CC11/TEMP/exampleunfilled_112RM_to_F99_NN_job.m';
            method={'interpolated',[]};
            fn=vertcat(fn{:});
            spm('Defaults','FMRI');
            addpath /imaging/dm01/fromRik/
            if ~exist(fullfile(groupdir,'PreOpt'),'dir'), mkdir(fullfile(groupdir,'PreOpt')); end
            if ~exist(fullfile(groupdir,'MaskROIs'),'dir'), mkdir(fullfile(groupdir,'MaskROIs')); end
            if ~exist(fullfile(groupdir,'PostOpt'),'dir'), mkdir(fullfile(groupdir,'PostOpt')); end
            if ~exist(fullfile(groupdir,'Mask-Pre'),'dir'), mkdir(fullfile(groupdir,'Mask-Pre')); end
            if ~exist(fullfile(groupdir,'Post-Mask'),'dir'), mkdir(fullfile(groupdir,'Post-Mask')); end
            if ~exist(fullfile(groupdir,'Post-Pre'),'dir'), mkdir(fullfile(groupdir,'Post-Pre')); end
            meanfile=cell(3,size(fn,1));
            for wit=1:3 % before/mask/after optimization
                fprintf('\nIteration stage %g',wit)
                tempfn=cell(1,size(fn{1,1},2));
                Vout=cell(1,size(fn{1,1},2));
                % average ROI connectivities per subject
                fprintf('\n\tAveraging ROI connectivity maps per subject')
                for sub=1:size(fn,1)
                    tempfnx=char(fn{sub,1}{wit,:});
                    spm_reslice(tempfnx,struct('mean',1,'which',0));
                    [tmppth, tmpnam]=fileparts(deblank(tempfnx(1,:)));
                    outmean=fullfile(tmppth,['mean' tmpnam '.nii']);
                    meanfile{wit,sub}=regexprep(outmean,'mean[^_]*_','meanAllROIs_');
                    movefile(outmean,meanfile{wit,sub});
                    fprintf('.')
                end % next sub
                for r=1:(size(fn{1,1},2) +1)
                    clear S
                    if r>size(fn{1,1},2)
                        tempfn{r}=char(meanfile(wit,:));
                    else
                        tempfn{r}=char(cellfun(@(x) x{wit,r},fn(:,1),'uniformoutput',0));
                    end
                    tempV=spm_vol(tempfn{r});
                    Y=spm_read_vols(tempV);
                    Vout{r}=rmfield(tempV(1),'pinfo');
                    [j, nam]=fileparts(Vout{r}.fname);
                    if wit==3, nam=regexprep(nam,'\d+$','Final'); end
                    nam=regexprep(nam,'_zC_','_C_');
                    Vout{r}.fname=fullfile(groupdir,[nam '.nii']);            
                    spm_write_vol(Vout{r},tanh(nanmean(Y,4)));
                    out=spm_jobman('serial', {RMSL_to_F99_job}, '',{Vout{r}.fname});
                    outmetricC(wit,r)=CaretScriptLinux_Nifti2Metric(out{1}.files{1},method{:}, true);
                    fprintf('.')
                    % rfx (this will obviously also calculate mean, so above is perhaps unnecessary)
                    % (but good to render C while doing stats on zC.)
                    switch wit
                        case 1, S.outdir = fullfile(groupdir,'PreOpt',nam);
                        case 2, S.outdir = fullfile(groupdir,'MaskROIs',nam);
                        case 3, S.outdir = fullfile(groupdir,'PostOpt',nam);
                    end
                    if ~exist(S.outdir,'dir'), mkdir(S.outdir); 
                    else delete(fullfile(S.outdir,'*.*'));
                    end
                    S.imgfiles{1}=cellstr(tempfn{r});
                    S.contrasts{1}.name = nam;
                    S.contrasts{1}.type = 'T';
                    S.contrasts{1}.c = 1; 
                     S.contrasts{2}.name=['Negative_' S.contrasts{1}.name];
                    S.contrasts{2}.type = 'T';
                    S.contrasts{2}.c = -1; 
                    batch_spm_anova(S);
                    threshoutname=localrenderspmcopy(S.outdir,aap);
                    for ton=1:numel(threshoutname)
                        out=spm_jobman('serial', {RMSL_to_F99_job}, '',threshoutname(ton));
                        outmetricCthresh(wit,r,ton)=CaretScriptLinux_Nifti2Metric(out{1}.files{1},method{:}, true);
                    end
                end % next roi
                %spm_reslice(char(cellfun(@(x) x.fname,Vout,'uniformoutput',0)),struct('mean',1,'which',0));
                fprintf('|')
            end % next iteration stage
            outgmean=cellstr(spm_select('FPList',groupdir,'^mean')); 
            for o=1:length(outgmean)
                try
                    movefile(outgmean{o},regexprep(outgmean{o},'mean[^_]*_','meanAllROIs_'))
                catch
                    % ignore this file if it already exists (will be overwritten) 
                end
            end
            
            %% now examine change from original > mask > optimised
            for wit={'Mask-Pre','Post-Mask','Post-Pre'}
                for r=1:size(fn{1,1},2) +1
                    switch wit{1};
                        case 'Mask-Pre', dims=[2 1]; ind=1;
                        case 'Post-Mask', dims=[3 2]; ind=2;           
                        case 'Post-Pre', dims=[3 1]; ind=3;               
                    end
                    if r>size(fn{1,1},2)
                        tempfn=meanfile(dims,:)';
                    else
                        tempfn=[cellfun(@(x) x{dims(1),r},fn(:,1),'uniformoutput',0),cellfun(@(x) x{dims(2),r},fn(:,1),'uniformoutput',0)];
                    end
                    for si=1:size(tempfn,1)
                        S.imgfiles{1}{si}=char(tempfn{si,:});
                    end
                    [j, nam]=fileparts(tempfn{1});
                    S.outdir=fullfile(groupdir,wit{1},nam);
                    if ~exist(S.outdir,'dir'), mkdir(S.outdir); 
                    else delete(fullfile(S.outdir,'*.*'));
                    end
                    S.contrasts{1}.name = [wit{1} '_' nam];
                    S.contrasts{1}.type = 'T';
                    S.contrasts{1}.c = [1 -1]; 
                    S.contrasts{2}.name=['Negative_' S.contrasts{1}.name];
                    S.contrasts{2}.type = 'T';
                    S.contrasts{2}.c = [-1 1]; 
                    batch_spm_anova(S);
                    threshoutname=localrenderspmcopy(S.outdir,aap)
                    for ton=1:numel(threshoutname)
                        out=spm_jobman('serial', {RMSL_to_F99_job}, '',threshoutname(ton));
                        outmetricCchangethresh(ind,r,ton)=CaretScriptLinux_Nifti2Metric(out{1}.files{1},method{:}, true);
                    end
                end
            end
            
            %% and iterated ROI overlap
            fprintf('\nCreating ROI overlap map')
            lastfiles=[lastfiles{:}]';
            clear allY
            for wit=[1 3] % before/after optimization
                for r=1:size(fn{1,1},2) % # seed ROIs
                    tempV=spm_vol(char(cellfun(@(x) x{r},lastfiles(:,wit),'uniformoutput',0)));
                    Y=spm_read_vols(tempV);
                    Vout=rmfield(tempV(1),'pinfo');
                    Vout.dt=[16 0];
                    [j, nam]=fileparts(Vout.fname);
                    if wit==3, nam=regexprep(nam,'\d+$','Final'); end
                    Vout.fname=fullfile(groupdir,[nam '.nii']);
                    allY(:,:,:,r,min(wit,2))=nanmean(Y,4); %(allY won't include masks, and values are proportion of subjects)
                    spm_write_vol(Vout,allY(:,:,:,r,min(wit,2)));
                    fprintf('.')
                    %%% now for Caret
                    out=spm_jobman('serial', {RMSL_to_F99_job}, '',{Vout.fname});
                    outmetric(wit,r)=CaretScriptLinux_Nifti2Metric(out{1}.files{1},method{:}, true);
                end
            end          
            allYa=nansum(allY(:,:,:,:,2),4);
            Vout.fname=fullfile(groupdir,[regexprep(nam,'.*_','allROIs_overlap') '.nii']);
            spm_write_vol(Vout,allYa);
            out=spm_jobman('serial', {RMSL_to_F99_job}, '',{Vout.fname});
            metriciteratedoverlap=CaretScriptLinux_Nifti2Metric(out{1}.files{1},method{:}, true);
            clear allYa
            threshY=uint16(allY);
            for r=1:size(fn{1,1},2)
                threshY(:,:,:,r,:)=threshY(:,:,:,r,:)*r;
            end
            threshY=nansum(threshY(:,:,:,:,2),4) + double(MPP{1}>0);
            Vout.fname=fullfile(groupdir,[regexprep(nam,'.*_','allROIs_ColouredAt50Percent') '.nii']);
            spm_write_vol(Vout,threshY);
            out=spm_jobman('serial', {RMSL_to_F99_job}, '',{Vout.fname});
            metricFinalROIsHalfSubs=CaretScriptLinux_Nifti2Metric(out{1}.files{1},'enclosing',[], true);
            threshY=allY*size(Y,4); threshY=uint8(round(threshY)==size(Y,4));
            for r=1:size(fn{1,1},2)
                threshY(:,:,:,r,:)=threshY(:,:,:,r,:)*r;
            end
            threshY=nansum(threshY(:,:,:,:,2),4) + double(MPP{1}>0);
            Vout.fname=fullfile(groupdir,[regexprep(nam,'.*_','allROIs_ColouredIfAllSubs') '.nii']);
            spm_write_vol(Vout,threshY);
            out=spm_jobman('serial', {RMSL_to_F99_job}, '',{Vout.fname});
            metricFinalROIsAllSubs=CaretScriptLinux_Nifti2Metric(out{1}.files{1},'enclosing',[], true);
            threshY=ceil(allY);
            for r=1:size(fn{1,1},2)
                threshY(:,:,:,r,:)=threshY(:,:,:,r,:)*r;
            end
            threshY=nansum(threshY(:,:,:,:,2),4) + double(MPP{1}>0);
            Vout.fname=fullfile(groupdir,[regexprep(nam,'.*_','allROIs_ColouredIfAnySubs') '.nii']);
            spm_write_vol(Vout,threshY);
            out=spm_jobman('serial', {RMSL_to_F99_job}, '',{Vout.fname});
            metricFinalROIsAnySubs=CaretScriptLinux_Nifti2Metric(out{1}.files{1},'enclosing',[], true);
            
            copyfile(fullfile(outdir{1},'MPP_iteration0.nii'),groupdir); 
            out=spm_jobman('serial', {RMSL_to_F99_job}, '',{fullfile(groupdir,'MPP_iteration0.nii')});
            metricInitialROIs=CaretScriptLinux_Nifti2Metric(out{1}.files{1},'enclosing',[], true);
            
            out=spm_jobman('serial', {RMSL_to_F99_job}, '',{fullfile(groupdir,'LabelledPriorRegions.nii')});
            metricROImasks=CaretScriptLinux_Nifti2Metric(out{1}.files{1},'enclosing',[], true);
            
            %%
            vv=spm_vol(char(regexprep(outmetricC(3,1:6),{'wF99','_mirrored.*'},{'','.nii'})));
            yy=spm_read_vols(vv);
            myy=min(yy,[],4);
            vvout=Vout;
            vvout.dt=[16 0];
            vvout.fname=fullfile(groupdir,'min_C_final.nii');
            spm_write_vol(vvout,myy)

            %% render for Caret
            addpath /imaging/dm01/AndrewBell/
           RenderCaretFiles(1,[],[3 2],struct(...
        'hemisphere','r','surf','flat',... %repmat({'flat'},[1 6])
        'metric',outmetricC(1,1:6),...
        'usealpha',0,...
        'view',[],... % repmat({''},[1 6])
        'caxis',[-0.3 0.3],...
        'func','',... % 'tanh(#)'
        'colormap',[],...
        'colorbar',0,'smooth',0,...
        'overlay','',...
        'overlayfunc','',... %round(#)/1
        'overlaycol',[0 1 0],...
        'title',''),...
        true) % final flag converts to bmp
        print('-dpng','-r82',fullfile(groupdir,'IteratedC_AllWindows_AllROIs_0p3.png'));

         RenderCaretFiles(1,[],[3 2],struct(...
        'hemisphere','r','surf','flat',... %repmat({'flat'},[1 6])
        'metric',outmetric(1,1:6),...
        'usealpha',1,...
        'view',[],... % repmat({''},[1 6])
        'caxis',[0 1],...
        'func','',... % 'tanh(#)'
        'colormap',[],...
        'colorbar',0,'smooth',0,...
        'overlay','',...
        'overlayfunc','',...
        'overlaycol',[0 1 0],...
        'title',''),...
        true) % final flag converts to bmp
        print('-dpng','-r82',fullfile(groupdir,'IteratedSeedOverlap_AllWindows_AllROIs_0p3.png'));
        
            %% also split ROIs by hemisphere for test data seed 2 seed matrix, and add any extra ROIs
            fprintf('\nCalculating seed-to-seed connectivity matrix for original ROI, optimised ROIs, and ROI masks:');
            tic;
            spmd
                
                S=struct([]);
                for type={'original','mask','optimised'}
                    S(1).(type{1}).seeds=roifilesO{labindex};
                    switch type{1}
                        case 'original'
                            S.(type{1}).seeds{1}=regexprep(S(1).(type{1}).seeds{1},'_window\d+_iteration\d+.nii','_iteration0.nii');
                        case 'optimised'
                            % leave as is
                        case 'mask'
                            S.(type{1}).seeds{1}=maskname;
                    end
                    aYlr=[];
                    labs={};
                    for rs=1:length(ROIs)
                        Yl=getdata(epinames{labindex},S(1).(type{1}).seeds{rs},'nanmean',[],GMmask{labindex},true,'l');
                        Yr=getdata(epinames{labindex},S(1).(type{1}).seeds{rs},'nanmean',[],GMmask{labindex},true,'r');
                        if labindex==1,fprintf('.'); end
                        Yl= R{labindex}{w}*Yl;
                        Yr= R{labindex}{w}*Yr;
                        aYlr=[aYlr Yl Yr];
                        labs=[labs strcat(S.(type{1}).seeds{rs},'_L') strcat(S.(type{1}).seeds{rs},'_R')];
                    end
                    %%%% seed to seed
                    Css=corr(aYlr);
                    Css=single(Css);
                    % corr can give nonsymmetric output due to floating point rounding errors!
                    Css=.5.*(Css+Css.');
                    Css=atanh(Css); % fisher transform
                    Css(isinf(Css))=max(Css(~isinf(Css)));
                    Css(1:(length(Css)+1):end)=0;
                    S.(type{1}).zCss(:,:)=Css;
                    [junk, temp]=spm_str_manip(labs,'C');
                    S.(type{1}).labs=temp.m;
                    %%%% can't save in spmd block
                    %save(fullfile(outdir{labindex},sprintf('Seed2Seed_window%g_%s.mat',w,type{1})),'zCss','seeds');
                end % next seed 'type'
                saveS(fullfile(outdir{labindex},sprintf('Seed2Seed_window%g.mat',w)),S);
            end % next subject
            fprintf(' (took %.2f mins)',toc/60)
            clear junk temp labs Css aYlr Yl Yr
            %             if strcmp(task,'doit')
            %                 figure(13); clf
            %                 positions=fitplots2(it);
            %                 try
            %                     for i=1:(it-1)
            %                         subplot('position',positions{i});
            %                         image(frame(i).cdata); axis off image;
            %                         title(sprintf('Iteration %g',i))
            %                     end
            %                 catch
            %                     % might not have bothered rendering?'
            %                 end
            %
            %                 print(10, fullfile(outdir,sprintf('IterationDiagnosticROIchange_window%g.png',w)),'-dpng');
            %                 print(11, fullfile(outdir,sprintf('IterationDiagnosticConnChange_window%g.png',w)),'-dpng');
            %             end
            
            %% plot connectivity matrices
            % (This is now done in PlotCmatForIteratedSSSVersion.m)
            sets={'original','mask','optimised'};
            nmain=length(maskname)*2;
            clear V D
            pos=fitplots2(4);
            for ss=1:length(S) % subjects
                temp=S{ss};
                for whichset=1:3 % iteration stages
                    temp2=temp.(sets{whichset}).zCss;
                    %%% hack to rearrange ROIs from alphabetical to John's favourite
                    reord=[3 2 5 4 1 6, 9 8 11 10 7 12, 13:(length(temp2))];
                    temp.(sets{whichset}).labs=temp.(sets{whichset}).labs(reord);
                    temp2=temp2(reord,reord);
                    %%%
                    M(1).(sets{whichset})(ss,:,:)=temp2; % full matrix
                    V(1).(sets{whichset})(ss,:)=squareform(temp2); % vectorised
                    temp2(1:(length(temp2)+1):end)=nan;
                    Rmain(ss,whichset)=nanmean(reshape(temp2(1:nmain,1:nmain),1,[])); % mean within MD
                    Rinter(ss,whichset)=nanmean(reshape(temp2(1:nmain,(nmain+1):end),1,[])); % mean within controls
                    Rcont(ss,whichset)=nanmean(reshape(temp2((nmain+1):end,(nmain+1):end),1,[])); % mean between MD and controls
                    D(ss,whichset)=Rmain(ss,whichset)-Rinter(ss,whichset); % difference of mean within MD and mean between 
                end
            end % next subject
            
            try delete(3); catch; end
            figure(3); clf; set(3,'color','w')
            cmax=-Inf;
            cmin=Inf;
            tax=nan(1,3);
            clear p
            for whichset=1:3
                %tax(whichset)=subplot(2,2,whichset);
                tax(whichset)=subplot('position',pos{whichset});
                toplot=squeeze(tanh(nanmean(M.(sets{whichset})))); % undo the fisher transform (after averaging across subjects)
                [h, p]=ttest(squeeze(V.(sets{whichset}))); % test each cell against zero
                h0=squareform(fdr_bh(p,0.05,'pdep'));
                if whichset>1
                    % for each cell, test increases and decreases from previous
                    [h, pinc]=ttest(squeeze(V.(sets{whichset})),squeeze(V.(sets{whichset-1})),0.025,'right');
                    [h, pdec]=ttest(squeeze(V.(sets{whichset})),squeeze(V.(sets{whichset-1})),0.025,'left');
                    hincdec=fdr_bh([pinc pdec],0.05,'pdep');
                    hinc=squareform(hincdec(1:length(h)));
                    hdec=squareform(hincdec((length(h)+1):end));
                    marks={{hinc,'+','color','r'},{hdec,'x','color','b'},{h0,'o','color','w'}};
                else
                    %%% just mark >0
                    %marks={h0,'o','color','w'};
                    %%% or include difference from previous (i.e. last)
                    [h, pinc]=ttest(squeeze(V.(sets{1})),squeeze(V.(sets{end})),0.025,'right');
                    [h, pdec]=ttest(squeeze(V.(sets{1})),squeeze(V.(sets{end})),0.025,'left');
                    hincdec=fdr_bh([pinc pdec],0.05,'pdep');
                    hinc=squareform(hincdec(1:length(h)));
                    hdec=squareform(hincdec((length(h)+1):end));
                    marks={{hinc,'+','color','r'},{hdec,'x','color','b'},{h0,'o','color','w'}};
                    
                end
                cmax=max(cmax,max(toplot(:)));
                cmin=min(cmin,min(toplot(:)));
                labs=temp.(sets{whichset}).labs;
                labs=regexprep(labs,{'.nii','.*/','_iteration\d+'},{'','',''});
                cla;
                try
                    PlotCMat(toplot,labs,'Correlation',marks);
                catch
                    imagesc(toplot); colorbar; axis image
                end
                title(sets{whichset})
            end
            for whichset=1:3
                %caxis(tax(whichset),([-cmax cmax]))
                %caxis(tax(whichset),([-floor(cmax*10)/10 cmax]))
                caxis(tax(whichset),([-0.35 0.6]))
                %caxis(tax(whichset),([cmin cmax]))
            end
            colormap jet
            colormap(flipud(lbmap(64,'redblue')))
            
            pos2=fitplots2([2 4]);
            tax(4)=subplot('position',pos2{6}); cla
            addpath /imaging/dm01/MoreTools/ % for PlotWithErrors
            opts.style='bars';
            opts.meanonly=1;
            opts.withinrows=-1;
            PlotWithErrors([Rmain(:,1), Rcont(:,1), Rinter(:,1)],['set','initial r',{'MD','Control','Inter'}],opts)
            txt=cell(1,3);
            [j, j, j, j, txt{1}]=ttestreport(Rmain(:,1));
            [j, j, j, j, txt{2}]=ttestreport(Rmain(:,1),Rcont(:,1));
            [j, j, j, j, txt{3}]=ttestreport(Rmain(:,1),Rinter(:,1));
            xlabel(txt);
            
            tax(5)=subplot('position',pos2{8}); cla
            opts.withinrows='previous';
            PlotWithErrors(D,['set','r(MD) - r(Inter)',sets],opts)
            txt=cell(1,3);
            [j, j, j, j, txt{1}]=ttestreport(D(:,1));
            [j, j, j, j, txt{2}]=ttestreport(D(:,1),D(:,2));
            [j, j, j, j, txt{3}]=ttestreport(D(:,2),D(:,3));
            xlabel(txt);
            
            fname=fullfile(groupdir,sprintf('Cmat_window%gof%g.eps',w,length(windowvec)));
            set(3,'renderer','painters');
            print(3,'-depsc',fname);
            %fix_pcolor_eps(fname); % very slow for large matrices

        end % next window
        
        diary off
        
        %%%% end here if just testing ROIs
        if ~strcmp(task,'doit'),
            resp=conn;
            return
        end
        %%%%
        
        fprintf('.Done.')
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end

return;

%%
function [TxS, allXYZvox, names, allXYZvox0]=getdata(imgs,rois,func,cov,mask,quiet,hem)
% inspired by JP_getdata on github; see also rex.m and spm_regions.m
tic

if ischar(rois), rois=cellstr(rois); end

if ~isstruct(imgs)
    if iscell(imgs), imgs=char(imgs); end
    imgs=spm_vol(imgs);
end

Vinv=inv(imgs(1).mat);

TxS=[]; allXYZvox=[]; names={}; allXYZvox0=[];

if exist('mask','var') && ~isempty(mask)
    Vm=spm_vol(mask);
    Ym=spm_read_vols(Vm);
    if size(Ym,4)>1
        Ym=prod(Ym,4);
    end
end

if ~exist('quiet','var') || isempty(quiet)
    quiet=0;
end

nvox=[];

for r=1:length(rois)
    
    V=spm_vol(rois{r});
    [pth, nam]=fileparts(V.fname);
    
    [allY, allXYZmm]=spm_read_vols(V);
    allY=round(allY);
    
    extratxt='';
    
    if exist('hem','var') && ischar(hem)
        switch lower(hem(1))
            case 'l'
                allY(allXYZmm(1,:)>0)=0;
                extratxt=[extratxt ' (left)'];
            case 'r'
                allY(allXYZmm(1,:)<0)=0;
                extratxt=[extratxt ' (right)'];
        end
    end
    
    if exist('mask','var') && ~isempty(mask)
        allY(~Ym)=0;
        extratxt=[extratxt ' (after masking)'];
    end
    
    vals=setdiff(unique(allY(:)),0);
    
    if isempty(vals)
        TxS=[TxS nan(length(imgs),1)];
        names=[names sprintf('%s_empty',nam)];
        if ~quiet
            fprintf('\n%s, 0 voxels! %s',nam,extratxt)
        end
    end
    
    for val=vals(:)'
        
        XYZmm=allXYZmm(:,allY==val);
        XYZvox0=Vinv(1:3,1:3)*XYZmm + repmat(Vinv(1:3,4),1,size(XYZmm,2));
        XYZvox=unique(round(XYZvox0'),'rows')';
        
        %%% if ROI is in different space to EPI, might need to fix rounding errors?
        trim=any((XYZvox-repmat(imgs(1).dim',1,size(XYZvox,2)))>0);
        XYZvox(:,trim)=[];
        %%%
        
        nvox(end+1)=size(XYZvox,2);
        if ~quiet
            fprintf('\n%s, val=%g: %g voxels %s',nam, val,nvox(end),extratxt)
        end
        
        Y=spm_get_data(imgs,XYZvox,0);
        
        if exist('func','var') && ~isempty(func)
            switch func
                case 'eigen' % code adapted from spm_regions
                    %                     % behzadi removed constant and linear trends per column before svd
                    %                     Y=Y-repmat(mean(Y),size(Y,1),1);
                    %                     Y=spm_detrend(Y,1);
                    %-Compute regional response in terms of first n eigenvariate
                    %--------------------------------------------------------------------------
                    [m n]   = size(Y);
                    if m > n
                        fprintf('\n%s (more TRs than voxels):',nam)
                        
                        % remove covariates?
                        if exist('cov','var') && ~isempty(cov)
                            R  = eye(size(Y,1)) - cov*pinv(cov);
                            Y = R*Y;
                            fprintf('Covariates removed;')
                        end
                        
                        YY=Y'*Y;
                        
                        [v s v] = svd(YY,0);
                        s       = diag(s); % eigenvalues
                        temp=s.^2; temp=cumsum(temp/sum(temp)); n=find(temp>=0.95,1);
                        v       = v(:,1:n); % eigenimage
                        u       = spm_en(Y*v); % eigenvariate
                    else
                        fprintf('\n%s (more voxels than TRs):',nam)
                        YY=Y*Y';
                        
                        % remove covariates? (borrowed from rex.m)
                        if exist('cov','var') && ~isempty(cov)
                            cov0=ones(size(YY,1),1);
                            cov1=detrend(cov,'constant');
                            proj=eye(size(YY,1))-[cov1,0*cov0]*pinv([cov1,cov0]); % removes covariates keeping scale unchanged
                            YY=proj*YY*proj';
                            fprintf('Covariates removed;')
                        end
                        
                        [u s u] = svd(YY,0);
                        s       = diag(s);
                        temp=s.^2; temp=cumsum(temp/sum(temp)); n=find(temp>=0.95,1);
                        u       = u(:,1:n); % eigenvariate
                        v       = spm_en(Y'*u); % eigenimage
                    end
                    fprintf('\t%g components explain %3g%% variance.',n,round(temp(n)*100));
                    d       = sign(sum(v));
                    Y       = u.*repmat(d,size(u,1),1);
                    Y=Y(:,1); fprintf(' Using first')
                otherwise
                    Y=eval(sprintf('%s(Y,2);',func));
                    
            end % how to summarise?
            
            % ???convert XYZvox to CoM here???
            
        end % summarise ROI?
        
        TxS=[TxS Y];
        allXYZvox=[allXYZvox XYZvox(:,:)];
        allXYZvox0=[allXYZvox0 XYZvox0(:,:)];
        names=[names sprintf('%s_%g',nam,val)];
        
    end % next roi label/value
    
    allXYZvox=single(allXYZvox);
    allXYZvox0=single(allXYZvox0);
end % next roi

if ~quiet
    fprintf(' (took %.2f mins)',toc/60);
    fprintf(' (Mean ROI size = %.1f (range %g-%g)',mean(nvox),min(nvox),max(nvox));
end

return

%%
function X=pinv_svds(A,varargin)
% djm: copied from pinv, but use svds instead of svd
[m,n] = size(A);

if n > m
    X = pinv(A',varargin{:})';
else
    %[U,S,V] = svd(A,0);
    [U,S,V] = svds(A,min(size(A))); % djm: the change!
    U(:,1)=-U(:,1); V(:,1)=-V(:,1); % djm: not sure this makes any difference
    if m > 1, s = diag(S);
    elseif m == 1, s = S(1);
    else s = 0;
    end
    if nargin == 2
        tol = varargin{1};
    else
        tol = max(m,n) * eps(max(s));
    end
    r = sum(s > tol);
    if (r == 0)
        X = zeros(size(A'),class(A));
    else
        s = diag(ones(r,1)./s(1:r));
        X = V(:,1:r)*s*U(:,1:r)';
    end
end
return

%%
function saveS(fname,S)
save(fname,'S');
return
