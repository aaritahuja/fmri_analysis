% AA module - write normalised EPIs
% Rhodri Cusack MRC CBU Cambridge Jan 2006-Aug 2007
% Resamples EPIs using *_seg_sn.mat file [if present] or *_sn.mat file
% Changed domain to once per session for improved performance when parallel

function [aap,resp]=aamod_norm_write(aap,task,i,j,currenttask)

if (ischar(aap))
    load(aap);
end;
if (ischar(i))
    i=str2num(i);
end;
if (ischar(j))
    j=str2num(j);
end;
if (exist('currenttask','var'))
    if (ischar(currenttask))
        tmp=load(currenttask);
        aap.tasklist.currenttask=tmp.currenttask;
    else
        aap.tasklist.currenttask=currenttask;
    end;
end;

resp='';

switch task
    case 'domain'
        resp='session';   % this module needs to be run once per subject

    case 'description'
        resp='SPM5 write normalised';

    case 'summary'
        sesspath=aas_getsesspath(i,j);
        resp=sprintf('Write normalised run on %s\n',sesspath);
    case 'report'
    case 'doit'

        tic
        
        % get the subdirectories in the main directory
        subj_dir = aas_getsubjpath(aap,i);

        subj.imgs = [];

        % Parameter file name
        structdir=fullfile(aas_getsubjpath(aap,i),aap.directory_conventions.structdirname);
        % changed to s* [de 210606]
        %subj.matname = dir(fullfile(structdir,['ms*' aap.acq_details.subjects(i).structuralfn '*seg_sn.mat']));
        % djm: just use suffix:
        subj.matname = spm_select('FPlist',structdir,[aap.acq_details.subjects(i).structuralfn '.*seg_sn.mat']);
        if isempty(subj.matname)
            subj.matname = spm_select('FPlist',structdir,'seg_sn.mat');
        end
        
        if (~length(subj.matname))
            subj.matname = dir(fullfile(structdir,['ms*' aap.acq_details.subjects(i).structuralfn '*_sn.mat']));
            if (~length(subj.matname))
                subj.matname = dir(fullfile(structdir,['s*' aap.acq_details.subjects(i).structuralfn '*_sn.mat']));
                if (~length(subj.matname))
                        aas_log(aap,1,sprintf('Cannot find normalisation _sn file in %s',fullfile(aas_getsubjpath(aap,i),aap.directory_conventions.structdirname)));
                end;
            end;
        end

        % Image to reslice
        P = aas_getimages(aap,i,j,aap.tasklist.currenttask.epiprefix);
        subj.imgs = strvcat(subj.imgs, P);
        % delete previous because otherwise nifti write routine doesn't
        % save disc space when you reslice to a coarser voxel
        for c=1:size(P,1)
            [pth fle ext]=fileparts(P(c,:));
            [s w]=aas_shell(['rm ' fullfile(pth,['w' fle ext])],true); % quietly
        end;
        
        
        % If we're doing the first session then...
        if (j==1)
            % Add t maps if they're being used
            P = dir(fullfile(subj_dir,aap.directory_conventions.tmapsdirname,'f*nii'));
            for b=1:length(P)
                subj.imgs = strvcat(subj.imgs, fullfile(subj_dir,aap.directory_conventions.tmapsdirname,P(b).name));
            end;
            % Add EPI mean?
            % djm: if this has not been realigned in same way as other epis, it's bounding box may
            % cause parts of other normalised images to be cropped. Safter not to include it here.
            %P=aas_getimages(aap,i,1,['mean']);
            %subj.imgs=strvcat(subj.imgs,P);
        end;
        % now write normalised
        if (length(subj.imgs)>0)
            % djm: just subj.matname
            spm_write_sn(subj.imgs,subj.matname,aap.spm.defaults.normalise.write);
        end;
        
        fprintf('\nTook %.1g minutes',toc/60)
        
    case 'checkrequirements'

    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;



