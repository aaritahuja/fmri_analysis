% AA module - segment and normalise structural, coregister EPI to this, and
% everything to template
% Daniel Mitchell 18/04/2012

% V2) Original version worked quite well, but BET not always successful and some GM possibly classified
% too much as CSF.
% Here I'll accept manually stripped brain if provided,
% (but this isn't accurate either, so try to coregister and normalise without it).
% Coregister EPIs to T2 if it exists rather than T1.
% Could also implement multimodal segmentation using T1 and T2, but this will take a bit of working
% out, and without EPI unsdistortion it's not clear it would be worth it.

% V2b)
% - reduce sampling for bias-correction step?
% - remove BET code?
% - Also segment T2 (if it exists), and average the tissue probabilities?
% - after initial coregistration and segmentation, use the segmentation to skull-strip both EPI and
%   structural, (then repeat coregistration of EPI to fine tune it)?
% - Illustrate coregistrations as tissue contours overlaid on EPI, T2, T1 and wT1?
% (ended up ignoring most of these except the new figure!)

% 2c) Try starting with the normalisation and segmentation, then coregister EPIs to native GM?
% Seems to work well. Perhaps the best I can do? inaccuracies seem due to EPI distortion
% (compare coregistered EPI to raw T1).

% 2d) Noticed that pre post lesion volumes don't line up well. Try to improve this using VBM8

function [aap,resp]=aamod_SegCoregNorm2d(aap,task,i)

resp='';

switch task
    case 'domain'
        resp='subject';   % this module needs to be run once per subject
    case 'whentorun'
        resp='justonce';  % should this be run everytime or justonce?
    case 'description'
        resp='Run initalprocessing';
    case 'summary'
        resp='coreg epi to structural; coreg all to template; normalise and segment';
    case 'report'
        
    case 'doit'
        
        subjpath=aas_getsubjpath(aap,i);
        switch aap.directory_conventions.T1template
            case '/imaging/dm01/AndrewBell/RhesusMacaqueAtlasTemplates/atlasfiles/112RM-SL_T1.nii'
                strucpath=spm_select('FPlist',fullfile(subjpath,'structurals'),'^mi\d+.*_structural.nii$');
            otherwise
                strucpath=spm_select('FPlist',fullfile(subjpath,'structurals'),'^s(CBU|MR).*.nii$');
        end
        [junk, T1name]=fileparts(strucpath);
        
        cwd=pwd;
        figure(spm_figure('GetWin'));
        
        %% Coregistration job configuration defaults
        % (normally just update mat without reslicing, but for T2, reslice to T1)
        
        clear coregjob
        coregjob{1}.spatial{1}.coreg{1}.estimate.other = {''};
        coregjob{1}.spatial{1}.coreg{1}.estimate.eoptions.cost_fun = 'nmi';
        switch aap.directory_conventions.T1template
            case '/imaging/dm01/AndrewBell/RhesusMacaqueAtlasTemplates/atlasfiles/112RM-SL_T1.nii'
                coregjob{1}.spatial{1}.coreg{1}.estimate.eoptions.sep = [4 2 1]; % default is [4 2]
                coregjob{1}.spatial{1}.coreg{1}.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001]./2; % half of the defaults
            otherwise
                coregjob{1}.spatial{1}.coreg{1}.estimate.eoptions.sep = [4 2]; % default is [4 2]
                coregjob{1}.spatial{1}.coreg{1}.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001]; 
        end
        coregjob{1}.spatial{1}.coreg{1}.estimate.eoptions.fwhm = [7 7];
        coregjob{1}.spatial{1}.coreg{1}.estimate.eoptions.graphics=0;
        
%         reslicejob{1}.spatial{1}.coreg{1}.estwrite.other = {''};
%         reslicejob{1}.spatial{1}.coreg{1}.estwrite.eoptions.cost_fun = 'nmi';
%         reslicejob{1}.spatial{1}.coreg{1}.estwrite.eoptions.sep = [4 2 1]; % default is [4 2]
%         reslicejob{1}.spatial{1}.coreg{1}.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001]./2; % half of the defaults
%         reslicejob{1}.spatial{1}.coreg{1}.estwrite.eoptions.fwhm = [7 7];
%         reslicejob{1}.spatial{1}.coreg{1}.estwrite.roptions.interp = 1;
%         reslicejob{1}.spatial{1}.coreg{1}.estwrite.roptions.wrap = [0 0 0];
%         reslicejob{1}.spatial{1}.coreg{1}.estwrite.roptions.mask = 1;
%         reslicejob{1}.spatial{1}.coreg{1}.estwrite.roptions.prefix = 'r'; % still prefixes r if this is set to empty
        
        %% normalise and segment settings
        defs = aap.spm.defaults.normalise;
        defs.estimate.weight = '';
        switch aap.directory_conventions.T1template
            case '/imaging/dm01/AndrewBell/RhesusMacaqueAtlasTemplates/atlasfiles/112RM-SL_T1.nii'
                defs.write.vox = [0.5 0.5 0.5];
                pth=fileparts(aap.directory_conventions.T1template);
                estopts.tpm=char({...
                    fullfile(pth,'gm_priors_ohsu+uw.nii,1'),...
                    fullfile(pth,'wm_priors_ohsu+uw.nii,1'),...
                    fullfile(pth,'csf_priors_ohsu+uw.nii,1'),...
                    });
                estopts.ngaus = [2 2 2 8]; % default = [2 2 2 4]
                estopts.biasfwhm= 30; % default = 60-75; reduce for monkey brain?
                estopts.samp = 2; % default=3; McLaren 2010 used 2
            otherwise
                defs.write.vox = [1 1 1];
                estopts.ngaus = [2 2 2 4]; 
                estopts.biasfwhm= 60; 
                estopts.samp = 3; 
        end
        defs.write.bb = aap.spm.defaults.normalise.write.bb;
        estopts.warpreg = 1; % default =  1
        estopts.warpco = 25; % default = 25; reduce later for monkey brain?
        % warpco can stay large here for speed
        estopts.biasreg = 0.0001; % default = 0.0001
        estopts.msk='';
        
        %% use normalise functions to create bias-corrected T1:
        % only write out attenuation corrected image
        writeopts.biascor = 1;
        writeopts.GM  = [0 0 0];
        writeopts.WM  = [0 0 0];
        writeopts.CSF = [0 0 0];
        writeopts.cleanup = 0;
        estopts.regtype=''; % turn off affine:
        %estopts.regtype='subj'; % turn on affine: subj corresponds to "average sized template" as used by McLaren 2010
        mT1name=fullfile(subjpath,'structurals',sprintf('m%s.nii',T1name));
        if ~exist(mT1name,'file'),
            fprintf('\nCreating bias-corrected T1...')
            warning off all
            out = spm_preproc(strucpath,estopts);
            warning on all
            [sn]   = spm_prep2sn(out); % don't write deformation files
            spm_preproc_write(sn,writeopts);
        else
            fprintf('\nBias-corrected T1 found.')
        end       
        
        %% bet: this is crude, but is just helping here with coregistration to template
        fprintf('\nSkull-stripping structural using BET...')
        sspath=strrep(mT1name,'.nii','_betT1.nii');
        cmd=sprintf('bet %s %s -v',mT1name,sspath);
        aas_runfslcommand(aap,cmd);
        % unzip everything
        cmd=sprintf('gunzip -rf %s',fileparts(strucpath));
        unix(cmd);
        
        %% coregister structurals to template, or to first structural of same individual
        % should perhaps have prededing dummy study module that forces this subject module to wait
        % for all structurals to be copied. In this analysis we know they've already been copied.
        try
            mname=regexprep(aap.acq_details.subjects(i).tag,'_.*','');
            otherm=strmatch(sprintf('%s_',mname),{aap.acq_details.subjects.tag})';
            if isempty(otherm), gocatch; end
        catch % no tags; coregister to template as normal
            otherm=i;
        end
        if i==otherm(1) % first scan from this individual, coregister to template
            target=aap.directory_conventions.T1template;
            fprintf('\nCoregistering to template')
        else
            % coregister to bias-corrected, coregistered first scan of this individual
            target=spm_select('FPList',fullfile(aas_getsubjpath(aap,otherm(1)),'structurals'),'^mm.*_structural.nii');
            if isempty(target)
                fprintf('\nWould like to coregister to first scan of this individual, which is still being processed.')
                fprintf('\nWill wait 3 minutes, then break and try again...')
                pause(180)
                error('...Deliberate error: waiting for first structural of this individual to be coregistered');
            else
                fprintf('\nCoregistering to first scan')
            end
        end
        
        ToTrans=[];
        switch aap.directory_conventions.T1template
            case '/imaging/dm01/AndrewBell/RhesusMacaqueAtlasTemplates/atlasfiles/112RM-SL_T1.nii'
                strucs=cellstr(spm_select('FPlist',fullfile(subjpath,'structurals'),...
                    sprintf('^m?%s.*.nii$',aap.directory_conventions.rawdataafterconversionprefix)));
            otherwise
                strucs=cellstr(spm_select('FPlist',fullfile(subjpath,'structurals'),...
                    sprintf('^m?s.*.nii$')));
        end
        ToTrans=[ToTrans; strucs];
        coregjob{1}.spatial{1}.coreg{1}.estimate.ref = {target};
        coregjob{1}.spatial{1}.coreg{1}.estimate.source = {sspath}; % {mT1name};
        coregjob{1}.spatial{1}.coreg{1}.estimate.other = ToTrans;
        fprintf('\nCoregistering structural images to template...')
        spm_jobman('run',coregjob)
        
        %% normalisation to McLaren template using attenuation corrected image
        switch aap.directory_conventions.T1template
            case '/imaging/dm01/AndrewBell/RhesusMacaqueAtlasTemplates/atlasfiles/112RM-SL_T1.nii'
                estopts.regtype='subj';    % turn on affine: subj corresponds to "average sized template" as used by McLaren 2010
                estopts.samp = 2; % default=3; McLaren 2010 used 2
                estopts.warpco = 12; % default = 25; reduce for monkey brain?
            otherwise
                estopts.regtype='mni';    % turn on affine: subj corresponds to "average sized template" as used by McLaren 2010
                estopts.samp = 3; % default=3; McLaren 2010 used 2
                estopts.warpco = 25; % default = 25; reduce for monkey brain?
        end
        estopts.warpreg = 1; % default =  1
        % (above two reduced from 1 to 0.1 and 12 to 5, to try better job of
        % 'correcting' lesion; didn't really work; could consider high-demensional warping?)
        % estopts.msk=spm_select('FPlist',fullfile(subjpath,'structurals'),'_mask');
        fprintf('\nNormalising & segmenting T1...')
        out = spm_preproc(mT1name,estopts);
        [sn,isn]=spm_prep2sn(out); % don't write deformation files yet...
        clear out
        
        writeopts.biascor = 0; % write native and modulated segmentations
        writeopts.GM  = [1 0 1];
        writeopts.WM  = [1 0 1];
        writeopts.CSF = [1 0 1];
        writeopts.cleanup = 1; % "light"
        spm_preproc_write(sn,writeopts);
        
        % now write forward and inverse deformation files
        snmatname=[spm_str_manip(mT1name,'sd') '_seg_sn.mat'];
        savefields(snmatname,sn);
        savefields([spm_str_manip(mT1name,'sd') '_seg_inv_sn.mat'],isn);
        
%         % (now try with vbm8 toolbox)
    % failed; perhaps because I don't have 6 tissue probability maps??
%          clear matlabbatch
%         matlabbatch{1}.spm.tools.vbm8.estwrite.data = {mT1name};
%         matlabbatch{1}.spm.tools.vbm8.estwrite.opts.tpm = estopts.tpm; % {'/imaging/local/software/spm_cbu_svn/releases/spm8_fil_r5236/toolbox/Seg/TPM.nii'};
%         matlabbatch{1}.spm.tools.vbm8.estwrite.opts.ngaus = [2 2 2 8]; % [2 2 2 3 4 2];
%         matlabbatch{1}.spm.tools.vbm8.estwrite.opts.biasreg = 0.0001; % might need to adjust this for clip artefacts, or try to mask them out?
%         matlabbatch{1}.spm.tools.vbm8.estwrite.opts.biasfwhm = 30; %60; % 
%         matlabbatch{1}.spm.tools.vbm8.estwrite.opts.affreg = 'subj'; % 'mni';
%         matlabbatch{1}.spm.tools.vbm8.estwrite.opts.warpreg = 1; % was 4; might need to reduce this to allow more severe warping?
%         matlabbatch{1}.spm.tools.vbm8.estwrite.opts.samp = 2; % could consider reducing this for accuracy?
%         matlabbatch{1}.spm.tools.vbm8.estwrite.extopts.dartelwarp.normlow = struct([]);
%         matlabbatch{1}.spm.tools.vbm8.estwrite.extopts.sanlm = 2;
%         matlabbatch{1}.spm.tools.vbm8.estwrite.extopts.mrf = 0.15;
%         matlabbatch{1}.spm.tools.vbm8.estwrite.extopts.cleanup = 1;
%         matlabbatch{1}.spm.tools.vbm8.estwrite.extopts.print = 0;
%         matlabbatch{1}.spm.tools.vbm8.estwrite.output.GM.native = 1;
%         matlabbatch{1}.spm.tools.vbm8.estwrite.output.GM.warped = 1;
%         matlabbatch{1}.spm.tools.vbm8.estwrite.output.GM.modulated = 2;
%         matlabbatch{1}.spm.tools.vbm8.estwrite.output.GM.dartel = 0;
%         matlabbatch{1}.spm.tools.vbm8.estwrite.output.WM.native = 1;
%         matlabbatch{1}.spm.tools.vbm8.estwrite.output.WM.warped = 1;
%         matlabbatch{1}.spm.tools.vbm8.estwrite.output.WM.modulated = 2;
%         matlabbatch{1}.spm.tools.vbm8.estwrite.output.WM.dartel = 0;
%         matlabbatch{1}.spm.tools.vbm8.estwrite.output.CSF.native = 1;
%         matlabbatch{1}.spm.tools.vbm8.estwrite.output.CSF.warped = 1;
%         matlabbatch{1}.spm.tools.vbm8.estwrite.output.CSF.modulated = 2;
%         matlabbatch{1}.spm.tools.vbm8.estwrite.output.CSF.dartel = 0;
%         matlabbatch{1}.spm.tools.vbm8.estwrite.output.bias.native = 1;
%         matlabbatch{1}.spm.tools.vbm8.estwrite.output.bias.warped = 1;
%         matlabbatch{1}.spm.tools.vbm8.estwrite.output.bias.affine = 0;
%         matlabbatch{1}.spm.tools.vbm8.estwrite.output.label.native = 1;
%         matlabbatch{1}.spm.tools.vbm8.estwrite.output.label.warped = 1;
%         matlabbatch{1}.spm.tools.vbm8.estwrite.output.label.dartel = 0;
%         matlabbatch{1}.spm.tools.vbm8.estwrite.output.jacobian.warped = 0;
%         matlabbatch{1}.spm.tools.vbm8.estwrite.output.warps = [0 0];
%         spm_jobman('run',matlabbatch)
        
        %% make maximum-liklihood tissue maps and brain mask (moved from aamod_makeROIs)
        fprintf('\nClassifying tissues...')
        gm=spm_select('FPlist',fullfile(subjpath,'structurals'),sprintf('^c1.*%s.*.nii$',T1name));
        wm=spm_select('FPlist',fullfile(subjpath,'structurals'),sprintf('^c2.*%s.*.nii$',T1name));
        csf=spm_select('FPlist',fullfile(subjpath,'structurals'),sprintf('^c3.*%s.*.nii$',T1name));
        Vg=spm_vol(gm); Yg=spm_read_vols(Vg); if isfield(Vg,'pinfo'),Vg=rmfield(Vg,'pinfo'); end
        Vw=spm_vol(wm); Yw=spm_read_vols(Vw); if isfield(Vw,'pinfo'),Vw=rmfield(Vw,'pinfo'); end
        Vc=spm_vol(csf); Yc=spm_read_vols(Vc);if isfield(Vc,'pinfo'),Vc=rmfield(Vc,'pinfo'); end
        Vg.fname=fullfile(aas_getsubjpath(aap,i),'structurals','GM.nii');
        spm_write_vol(Vg,Yg>Yw & Yg>Yc & Yg>0.05);
        Vw.fname=fullfile(aas_getsubjpath(aap,i),'structurals','WM.nii');
        spm_write_vol(Vw,Yw>Yg & Yw>Yc & Yw>0.05);
        Vc.fname=fullfile(aas_getsubjpath(aap,i),'structurals','CSF.nii');
        spm_write_vol(Vc,Yc>Yw & Yc>Yg & Yc>0.05);
        clear Yg Yw Yc
        
        % make brain mask
        segmentations=spm_select('FPlist',fullfile(subjpath,'structurals'),sprintf('^c\\d.*%s.*.nii$',T1name));
        V=spm_vol(segmentations);
        Y=spm_read_vols(V);
        if isfield(V,'pinfo'),V=rmfield(V,'pinfo'); end
        V(1).fname=fullfile(subjpath,'structurals','ssbrainmask.nii');
        if isfield(V,'pinfo'); V=rmfield(V,'pinfo'); end
        brainmask=sum(Y,4)>0.5;
        brainmask=spm_dilate(double(brainmask));
        spm_write_vol(V(1),brainmask);
        
        %% and write normalised T1 and segmentations
        spm_write_sn(char(segmentations,mT1name,Vg.fname,Vw.fname,Vc.fname,V.fname),snmatname,defs.write);
        
        %% skull-strip T1 (nicer viewing without hyperintensity at ears)
        V=spm_vol(mT1name); Y=spm_read_vols(V);
        if isfield(V,'pinfo'),V=rmfield(V,'pinfo'); end
        V.fname=fullfile(subjpath,'structurals','ssT1.nii');
        spm_write_vol(V,Y.*brainmask);
        clear Y brainmask
        
        %% coregister epis to native GM segmentation
        %%
        % (The reverse is potentially a bit quicker, but I'd like to get everything coregistered
        % with template, and there's no need to reslice at the moment. T2 should be more accurate
        % as it has similar contrast. Non-bias-corrected should be more accurate as EPIs are not
        % bias-corrected?)
        ToTrans=cellstr(spm_select('FPlist',aas_getsesspath(aap,i,1),sprintf('^%s%s.*.nii$',aap.tasklist.currenttask.epiprefix,aap.directory_conventions.rawdataafterconversionprefix)));
        anEPI=ToTrans{round(length(ToTrans)/2)};
        
        % IN CASE WANT TO TEST DIFFERENT OPTIONS:
        % anEPI=regexprep(anEPI,'.nii','_temp.nii'); copyfile(ToTrans{end},anEPI); ToTrans=anEPI;
        coregjob{1}.spatial{1}.coreg{1}.estimate.ref = {gm};
        coregjob{1}.spatial{1}.coreg{1}.estimate.source = {anEPI};
        coregjob{1}.spatial{1}.coreg{1}.estimate.other = ToTrans;
        fprintf('\nCoregistering EPIs to structural...')
        spm_jobman('run',coregjob)
        
        %% Now save graphical check
        fprintf('\nPlotting diagnostic...')
        H=spm_figure('GetWin');
        figure(H); clf; set(H,'color','k');
        % use ss structurals to avoid intensity hotspots on ears which make brain relatively dark
        sliceperview=[10, 2, 0]; % try to avoid midline fissure
            underlays={anEPI,...
                fullfile(subjpath,'structurals','ssT1.nii'),...
                strrep(mT1name,'structurals/','structurals/w'),...
                aap.directory_conventions.T1template};

        for base=1:length(underlays)
            prop=1.5;
            if base<(length(underlays)-1)
                %so=slover(char(underlays{base},Vg.fname,Vw.fname,Vc.fname));%
                so=slover(char(Vg.fname,underlays{base}));%
                sc=1.7;
                views={spm_matrix([0 0 0 0 0 0 sc sc sc]),...
                    spm_matrix([-25 0 0 pi/2 0 -pi/2 -sc -sc sc]),...
                    spm_matrix([0 0 0 pi/2 0 0 sc -sc sc])};
            else
                %                 so=slover(char(underlays{base}, ...
                %                     strrep(Vg.fname,'structurals\','structurals\w'), ...
                %                     strrep(Vw.fname,'structurals\','structurals\w'), ...
                %                     strrep(Vc.fname,'structurals\','structurals\w')));%
                so=slover(char(strrep(Vg.fname,'structurals/','structurals/w'),underlays{base}));%
                sc=0.9;
                views={spm_matrix([0 0 0 0 0 0 sc sc sc]),...
                    spm_matrix([-12 0 0 pi/2 0 -pi/2 -sc -sc sc]),...
                    spm_matrix([0 0 0 pi/2 0 0 sc -sc sc])};
                if base==length(underlays), prop=0.9; end
            end
            
            so.img(2).prop=prop;so.img(2).cmap=colormap('gray'); %so.img(1).func='mmax=prctile(i1(:),95); i1(i1>mmax)=mmax;i1=i1./mmax*100;'; %so.img(1).func='i1=i1./median(i1(:))*100;';so.img(1).range=[0 200];
            %so.img(2).type='contour'; so.img(2).contours=[eps eps]; so.img(2).linespec='g-'; so.img(2).linewidth=1;
            %so.img(3).type='contour'; so.img(3).contours=[eps eps]; so.img(3).linespec='r-'; so.img(3).linewidth=1;
            %so.img(4).type='contour'; so.img(4).contours=[eps eps]; so.img(4).linespec='c-'; so.img(4).linewidth=1;
            %so.img(2).prop=0.1;so.img(2).cmap=colormap('autumn'); so.img(2).range=[0.5 1]; so.img(2).outofrange={[],[]};
            so.img(1).prop=0.2;so.img(1).cmap=repmat([1 0 0],[128 1]); so.img(1).range=[0.5 2]; so.img(1).outofrange={[],[]};
            so.figure=spm_figure('GetWin');
            so.cbar=[];
            so.clf=0;
            so.labels='none';
            for col=1:3
                so.transform=views{col};
                so.area.position=[(col-1)/3 (base-1)/5 1/3 1/5];
                so.slices=sliceperview(col);%5:5:30;%
                so=paint(so);
            end
            annotation('textbox',so.area.position,'string',regexprep(so.img(2).vol.fname,'.*/',''),...
                'LineStyle','none','color','w','horizontalAlignment','right');
        end
        annotation('textbox',[0 0.95 1 0.05],'string',subjpath,'LineStyle','none','color','r');
        print('-dpng',fullfile(subjpath,'diagnostic_aamod_SegCoregNorm2d.png'));
        set(H,'color','w');
        
        %% reslice native structurals?
        %spm_reslice(char(aap.directory_conventions.T1template,char(strucs)),struct('mean',0,'which',1,'prefix','r'))
        
        cd(cwd);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;

%------------------------------------------------------------------------
function savefields(fnam,p)
if length(p)>1, error('Can''t save fields.'); end;
fn = fieldnames(p);
if numel(fn)==0, return; end;
for i=1:length(fn),
    eval([fn{i} '= p.' fn{i} ';']);
end;
if str2double(version('-release'))>=14,
    save(fnam,'-V6',fn{:});
else
    save(fnam,fn{:});
end;

return;
%------------------------------------------------------------------------