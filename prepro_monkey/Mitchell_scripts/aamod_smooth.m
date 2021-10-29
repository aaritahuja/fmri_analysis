% AA module - smoothing
% Kernel size determined by aap.spm_analysis.FWHM
% Rhodri Cusack MRC CBU Cambridge Jan 2006- Aug 2007
% Now once per session for improved performance when running in parallel
%
% djm: added option to restrict smoothing to be within masks (specified as cell array of files)
% see e.g.
% https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind1008&L=spm&P=R45120&1=spm&9=A&J=on&K=2&d=No+Match%3BMatch%3BMatches&z=4
% https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind0901&L=spm&P=R14635&1=spm&9=A&J=on&K=3&d=No+Match%3BMatch%3BMatches&z=4

function [aap,resp]=aamod_smooth(aap,task,i,j,currenttask)

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
        resp='session';   % this module needs to be run once per session
        
    case 'description'
        resp='SPM5 smooth';
        
    case 'summary'
        subjpath=aas_getsubjpath(i);
        resp=sprintf('Smooth run on %s\n',subjpath);
    case 'report'    
    
    case 'doit'
        
        % get the subdirectories in the main directory
        subj_dir = aas_getsubjpath(aap,i);
        
        subj.imgs = [];
        
        % Add t maps if they're being used
        P = dir(fullfile(subj_dir,'functionals',aap.directory_conventions.tmapsdirname,'w*nii'));
        for b=1:length(P)
            subj.imgs = strvcat(subj.imgs, fullfile(subj_dir,aap.directory_conventions.tmapsdirname,P(b).name));
        end;
        % Image to reslice
            % get files in this directory
            P = aas_getimages(aap,i,j,[aap.tasklist.currenttask.epiprefix]);
            subj.imgs = strvcat(subj.imgs, P);
                
        % now smooth    
        s   = aap.tasklist.currenttask.settings.FWHM;
        n   = size(subj.imgs,1);
        
        masks=aap.tasklist.currenttask.settings.masks;
        if ~isempty(masks)
            fprintf('\nSmoothing masks')
            Vm=cell(1,length(masks));
            Ym=cell(1,length(masks));
            Yms=cell(1,length(masks));
            reslicedmask=cell(1,length(masks));
            for m=1:length(masks);
                anEPI=deblank(subj.imgs(1,:));
                
                [pth nam]=fileparts(masks{m});
                thismask=spm_select('FPList',fullfile(subj_dir,pth),nam);
                spm_reslice(char(anEPI,thismask),struct('which',1,'mean',0))
                reslicedmask{m}=regexprep(thismask,'(.*)/(.*)','$1/r$2');
                
                Vm{m}=spm_vol(reslicedmask{m});          
                Ym{m}=round(spm_read_vols(Vm{m}));
                Yms{m}=nan(size(Ym{m}));
                [BB, vx{m}] = spm_get_bbox(Vm{m});
                spm_smooth(Ym{m}, Yms{m},abs(s./[vx{m}(1) vx{m}(2) vx{m}(3)]));
                Ym{m}=logical(Ym{m});
                fprintf('.')
            end
        end
        
        spm_progress_bar('Init',n,'Smoothing','Volumes Complete');
        fprintf('\nSmoothing EPIs')
        warning off all % don't worry about divide by zero
        for i = 1:n
            Q = deblank(subj.imgs(i,:));
            [pth,nam,xt,nm] = spm_fileparts(Q);
            U = fullfile(pth,[sprintf('s%g',s) nam xt nm]); %djm: add smoothing size
            if isempty(masks)
                spm_smooth(Q,U,s);
            else % djm: smooth within masks, but not between them
                VQ=spm_vol(Q);
                YQ=spm_read_vols(VQ);
                YU=YQ;
                for m=1:length(masks);
                    YQs=nan(size(YQ));
                    spm_smooth(YQ.*Ym{m},YQs,abs(s./[vx{m}(1) vx{m}(2) vx{m}(3)]));
                    YQs=YQs./Yms{m}; % this 'undoes' the zero padding of convolution at mask edges
                    YU(Ym{m})=YQs(Ym{m});
                end
                VQ.fname=U;
                spm_write_vol(VQ,YU);
            end
            spm_progress_bar('Set',i);
            if ~mod(i,50), fprintf('.'); end
        end
        spm_progress_bar('Clear');
        warning on all 
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;



