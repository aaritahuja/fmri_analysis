% AA module
%
function [aap,resp]=aamod_makeROIs(aap,task,i)

resp='';

switch task
    case 'domain'
        resp='subject';   % this module needs to be run once per subject

    case 'description'
        resp='customised script to make ROIs';

    case 'summary'
        subjpath=aas_getsubjpath(aap,i);
        resp=sprintf('ROIs for %s\n',subjpath);

    case 'report'

    case 'doit'

        %imcalcflags={0,0,16,0};

        templaterois=aap.tasklist.currenttask.settings.templaterois;

        if exist(templaterois,'dir')
            templaterois=cellstr(spm_select('FPlist',templaterois,'.*.nii$'));
        elseif ischar(templaterois)
            templaterois=cellstr(templaterois);
        end
        
        %%%%% create subject roi dir if necessary
        roidir=fullfile(aas_getsubjpath(aap,i),'rois');
        if ~exist(roidir,'dir'), mkdir(roidir); end
        structdir=fullfile(aas_getsubjpath(aap,i),aap.directory_conventions.structdirname);

        %%% make range of possible tissue masks
        fprintf('\nMaking tissue masks')
        tissuenames={'c1m','c2m','c3m'};

        if isfield(aap.tasklist.currenttask,'epiprefix') && ~isempty(aap.tasklist.currenttask.epiprefix)
           epifile=aas_getimages(aap,i,1,aap.tasklist.currenttask.epiprefix,aap.acq_details.numdummies+1);
           if isempty(strfind(aap.tasklist.currenttask.epiprefix,'w'))
               % reslice the unnormalised strutcurals (because the epi has not been normalised; assume coregistered)
               tis=spm_select('FPlist',structdir,'^c.*');
               normprefixes={'','w','r'};
           else
               % reslice the normalised strutcurals
               tis=spm_select('FPlist',structdir,'^wc.*');
               normprefixes={'','w','rw'};
           end
           spm_reslice(char(epifile,tis),struct('which',1,'mean',0)); % linear interpolation
           
        else
            normprefixes={'','w'};
        end
        
        for normprefix=normprefixes
            strucfile=spm_select('FPlist',structdir,sprintf('^%smm.*structural.nii',normprefix{1}));
            if numel(strucfile)==1 || isempty(strucfile)
                strucfile=spm_select('FPlist',structdir,sprintf('^%ss.*structural.nii',normprefix{1}));
            end
            if numel(strucfile)==1 || isempty(strucfile)
                strucfile=spm_select('FPlist',structdir,sprintf('^%sms.*structural.nii',normprefix{1}));
            end
            if numel(strucfile)==1 || isempty(strucfile)
                strucfile=spm_select('FPlist',structdir,sprintf('^%smm.*structural.nii',normprefix{1}(2:end)));
            end
            if numel(strucfile)==1 || isempty(strucfile)
                strucfile=spm_select('FPlist',structdir,sprintf('^%sms.*.nii',regexprep(normprefix{1},'^r',''))); % for human?
                strucfile=deblank(strucfile(1,:));
            end

            masks={};
            masks_processed={};
            for class=tissuenames
                V=spm_vol(spm_select('FPlist',structdir,sprintf('^%s%s.*.nii',normprefix{1},class{1})));
                Y=spm_read_vols(V);
                origname=V.fname;
                
                for t=[50 90 99];
                    %V.fname=regexprep(origname,'c(\d)m',[sprintf('t%g',t),'c$1m']);
                    V.fname=regexprep(origname,'(.*)/(.*)',['$1/' sprintf('t%g',t),'$2']);
                    masks{end+1}=V.fname;
                    T=Y>(t/100);
                    spm_write_vol(V,T);
                    %masks{end+1}=V.fname;
                    fprintf('.')
                    % may want to tidy up here by growing/eroding/cleaning; also select processed
                    % masks to show
                    if strcmp(class{1},'c1m'),
                        if t==90, masks_processed{end+1}=V.fname; end
                    elseif strcmp(class{1},'c2m'),
                        % double erosion as in Behzadi (probably too severe for small monkey brain?)
                        T=spm_erode(double(T));
                        V.fname=regexprep(V.fname,'(t\d+)','e$1');
                        spm_write_vol(V,T);
                        if t==99, masks_processed{end+1}=V.fname; end
                         T=spm_erode(double(T));
                        V.fname=regexprep(V.fname,'(et\d+)','e$1');
                        spm_write_vol(V,T);                      
                    elseif strcmp(class{1},'c3m'),
                        % remove singleton voxels (not clear exactly what behzadi used)   
                        T=bwareaopen(T,2);
                        V.fname=regexprep(V.fname,'(t\d+)','o$1');
                        spm_write_vol(V,T);
                        if t==99, masks_processed{end+1}=V.fname; end
                    end
                    fprintf('.')
                    end
            end
            % diagnostic (all thresholds, no morphological processing)
            so=slover(char(strucfile,char(masks)));
            so.img(1).prop=1;
            for r=1:length(so.img)
                if ~isempty(strfind(so.img(r).vol.fname,'c1m'))
                    so.img(r).cmap=repmat([0 1 0],64,1); so.img(r).range=[0.5 1]; so.img(r).outofrange={0,1}; so.img(r).prop=0.15;
                elseif ~isempty(strfind(so.img(r).vol.fname,'c2m'))
                    so.img(r).cmap=repmat([1 0 0],64,1); so.img(r).range=[0.5 1]; so.img(r).outofrange={0,1}; so.img(r).prop=0.15;
                elseif ~isempty(strfind(so.img(r).vol.fname,'c3m'))
                    so.img(r).cmap=repmat([0 0.6 1],64,1); so.img(r).range=[0.5 1]; so.img(r).outofrange={0,1}; so.img(r).prop=0.25;
                else
                    so.img(r).range=so.img(r).range./2;
                end
            end
            so.figure=spm_figure('GetWin'); so.clf=1;
            so.area.position=[0 0 1 1];
            so.slices=5:5:30;
            so.cbar=[];
            %so.labels.format='%3.0fmm';
            so=paint(so);
            print(so.figure,'-dpng',fullfile(aas_getsubjpath(aap,i),sprintf('diagnostic_TissueMasks_%s.png',normprefix{1})));
            
            % repeat highest threshold, after morphological processing
            so=slover(char(strucfile,char(masks_processed)));
            so.img(1).prop=1;
            for r=1:length(so.img)
                if ~isempty(strfind(so.img(r).vol.fname,'c1m'))
                    so.img(r).cmap=repmat([0 1 0],64,1); so.img(r).range=[0.5 1]; so.img(r).outofrange={0,1}; so.img(r).prop=0.4;
                elseif ~isempty(strfind(so.img(r).vol.fname,'c2m'))
                    so.img(r).cmap=repmat([1 0 0],64,1); so.img(r).range=[0.5 1]; so.img(r).outofrange={0,1}; so.img(r).prop=0.4;
                elseif ~isempty(strfind(so.img(r).vol.fname,'c3m'))
                    so.img(r).cmap=repmat([0 0.6 1],64,1); so.img(r).range=[0.5 1]; so.img(r).outofrange={0,1}; so.img(r).prop=0.5;
                else
                    so.img(r).range=so.img(r).range./2;
                end
            end
            so.figure=spm_figure('GetWin'); so.clf=1;
            so.area.position=[0 0 1 1];
            so.slices=5:5:30;
            so.cbar=[];
            %so.labels.format='%3.0fmm';
            so=paint(so);
            print(so.figure,'-dpng',fullfile(aas_getsubjpath(aap,i),sprintf('diagnostic_TissueMasks_%s_processed.png',normprefix{1})));
        end

        %%%%% back-normalise template-space ROIs
        if isempty(templaterois)
            fprintf('\nNo template ROIs provided.')
        else
            fprintf('\nCreating & back normalising template ROIs')

            % Images to unnormalise (copied to subject roi dir)
            clear P
            warning off all; rmpath /imaging/local/spm/spm5/toolbox/Anatomy; warning on all

            % copy template ROIs to local folder
            localrois=cell(size(templaterois));
            for r=1:length(templaterois)
                [opth onam oext]=fileparts(templaterois{r});
                localrois{r}=fullfile(roidir,[onam oext]);
                copyfile(templaterois{r},localrois{r});              
                if ~mod(r,2), fprintf('.'); end
            end

            cd(roidir);

            %%% reslice template ROIs to space of normalised brain
            %filt=['^' regexprep(aap.tasksettings.aamod_norm_write.outputstreams.stream.filter,'\*','m.*') '$'];
            filt='^wm.*nii$';
            strucfile=spm_select('FPlist',structdir,filt);
            resliceflags = struct('interp',1,'mask',0,'mean',0,'which',1,'wrap',[0 0 0]');
            spm_reslice(char(strucfile,char(localrois)),resliceflags)

            if aap.tasklist.currenttask.settings.warptonative
                %%% write unnormalised
                invmatname = spm_select('FPlist',structdir,'.*seg_inv_sn.mat');
                %flags=struct('bb',spm_get_bbox(aap.directory_conventions.T1template),'vox',[0.5 0.5 0.5]);
                flags=struct('bb',spm_get_bbox(spm_select('FPlist',structdir,'^mmi.*structural.nii$')),'vox',[0.5 0.5 0.5]);
                spm_write_sn(char(strcat('r',templaterois)),invmatname,flags); %,defs.write
            end
            
            % conn may want rois to match structural; so far this is true in normalised space, but
            % would need to be fixed for native space.
            %spm_check_orientations(spm_vol(char(spm_select('FPlist',structdir,'^w?mm.*.nii'),spm_select('FPlist',roidir,sprintf('^(w|r)%s',templaterois{1})))))

            % also mask all ROIs to include only >x% GM
            % note that physiological aftefacts elsewhere are only a problem if they also exist in the
            % seed, so perhaps good to be conservative here? But might be
            % better/easier to do the musking in the connectivity script. 
            th=90;
            fprintf('\nMasking ROIs to GM:')
            allROIs=cellstr(spm_select('list',roidir,'.*.nii$'));
            nGM=spm_select('FPlist',structdir,sprintf('^t%gc1.*.nii$',th));
            wGM=spm_select('FPlist',structdir,sprintf('^t%gwc1.*.nii$',th));
            YnGM=spm_read_vols(spm_vol(nGM));
            YwGM=spm_read_vols(spm_vol(wGM));
            if aap.tasklist.currenttask.settings.warptonative
                prefixes={'wr','r'};
            else
                prefixes={'r'};
            end
            for normprefix=prefixes; %'wr' w=back-normalized to native space
                labs={};
                for r=1:length(localrois)
                    V=spm_vol(regexprep(localrois{r},'(.*)/',sprintf('$1/%s',normprefix{1})));
                    Y=spm_read_vols(V);
                    %V.fname=fullfile(roidir,[sprintf('GM%g',th) normprefix{1}, templaterois{r}]);
                    V.fname=regexprep(localrois{r},'(.*)/',sprintf('$1/GM%g%s',th,normprefix{1}));
                    if strcmp('spm - 3D normalized',V.descrip)
                        M=round(Y.*YnGM);
                    else
                        M=round(Y.*YwGM);
                    end

                    %fprintf('\n%s: %g voxels',V.fname,sum(M(:)))
                    U=unique(M(:));
                    for u=1:length(U);
                        if sum(M(:)==U(u))>0,
                            fprintf('\n%s: %g voxels of value %g',V.fname,sum(M(:)==U(u)),U(u))
                            labs{r}=sprintf('\n%s:\t%g vox of value %g',regexprep(V.fname,{'.*/','\..*'},{'',''}),sum(M(:)==U(u)),U(u));
                        else
                            M(Y(:)==U(u))=U(u);
                            fprintf('\n%s: NO VOXELS OF VALUE %G SURVIVE MASK, USING ALL %g UNMASKED VOXELS',V.fname,U(u),sum(Y(:)==U(u)))
                            labs{r}=sprintf('\n%s:\tNO VOXELS OF VALUE %G SURVIVE MASK, USING ALL %g UNMASKED VOXELS',regexprep(V.fname,{'.*/','\..*'},{'',''}),U(u),sum(Y(:)==U(u)));
                        end
                    end
                    spm_write_vol(V,M);
                end
            end

            %%%%% render and print diagnostic for GM-masked ROIs in template space
            fprintf('\nPrinting diagnostic:')
            % diagnostic
            cols=hsv(length(templaterois));
            strucfile=spm_select('FPlist',structdir,sprintf('^wmm.*.nii'));
            if isempty(strucfile), 
                strucfile=spm_select('FPlist',structdir,sprintf('^wms.*.nii')); % human?
            end
            so=slover(char(strucfile,char(regexprep(localrois,'(.*)/',sprintf('$1/GM%gr',th)))));
            so.img(1).prop=1; so.img(1).range=so.img(1).range./2;
            for r=2:length(so.img)
                so.img(r).cmap=repmat(cols(r-1,:),64,1);
                so.img(r).range=[0.5 1];
                so.img(r).outofrange={0,1};
                so.img(r).prop=1;
            end
            so.figure=spm_figure('GetWin'); so.clf=1;
            so.area.position=[0 0 1 1];
            so.slices=5:2:35;
            so.cbar=[];
            so=paint(so);
            H=[];
            if length(templaterois)<50
                for r=1:length(templaterois)
                    H(end+1)=annotation('textbox',[0 (r)*0.02 1 0.02],'string',labs{r},'color',mean([.75 .75 .75;cols(r,:)]),'edgecolor','none');
                end
            end
            print(so.figure,'-dpng',fullfile(aas_getsubjpath(aap,i),sprintf('diagnostic_wGM%g_ROIs.png',th)));
            delete(H)
        end

    case 'checkrequirements'

    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;

return


