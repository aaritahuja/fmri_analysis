% AA module - realignment
% Realignment using SPM5
% i=subject num
% Rhodri Cusack MRC CBU 2004-6 based on original by Matthew Brett
% Modified by Rik Henson 2006-8 to accept reslice "which" option 
% 	(plus more defaults can be passed)

function [aap,resp]=aamod_realign(aap,task,i,currenttask)
if (ischar(aap))
    load(aap);
end;
if (ischar(i))
    i=str2num(i);
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
        resp='subject';
    case 'description'
        resp='SPM5 realignment';
    case 'summary'
        resp='Done SPM5 realignment\n';
    case 'report'
        mvmean=[];
        mvmax=[];
        mvstd=[];
        mvall=[];
        nsess=length(aap.acq_details.sessions);

        qq=[];
        for j=1:nsess
            imfn=aas_getimages(aap,i,j,aap.tasklist.currenttask.epiprefix);
            imV=spm_vol(imfn);
            mv=[];
            for k=1:length(imV)
                tmp=spm_imatrix(imV(k).mat);
                mv=[mv;tmp(1:6)];
            end;
            if (j==1)
                firstscan=mv(1,:);
            end;
            mv=mv-repmat(firstscan,[size(mv,1) 1]);

            mv(:,4:6)=mv(:,4:6)*180/pi; % convert to degrees!
            mvmean(j,:)=mean(mv);
            mvmax(j,:)=max(mv);
            mvstd(j,:)=std(mv);
            mvall=[mvall;mv];
        end;

        
        aap.report.html=strcat(aap.report.html,'<h3>Movement maximums</h3>');
        aap.report.html=strcat(aap.report.html,'<table cellspacing="10">');
        aap.report.html=strcat(aap.report.html,sprintf('<tr><td align="right">Sess</td><td align="right">x</td><td align="right">y</td><td align="right">z</td><td align="right">rotx</td><td align="right">roty</td><td align="right">rotz</td></tr>',j));
        for j=1:nsess
        aap.report.html=strcat(aap.report.html,sprintf('<tr><td align="right">%d</td>',j));
        aap.report.html=strcat(aap.report.html,sprintf('<td align="right">%8.3f</td>',mvmax(j,:)));
        aap.report.html=strcat(aap.report.html,sprintf('</tr>',j));
        end;
        aap.report.html=strcat(aap.report.html,'</table>');
        
        
        
        
        
        aap=aas_report_addimage(aap,fullfile(aas_getsubjpath(aap,i),'diagnostic_aamod_realign.jpg'));     

    case 'doit'
        % Get realignment defaults
        defs = aap.spm.defaults.realign;
        
        % JC SPM 8 update fix 190411
        if ~isfield(defs.estimate,'weight')
            defs.estimate.weight = [];
        end
                
        % Flags to pass to routine to calculate realignment parameters
        % (spm_realign)
        reaFlags = struct(...
            'quality', defs.estimate.quality,...  % estimation quality
            'fwhm', defs.estimate.fwhm,...        % smooth before calculation
            'rtm', defs.estimate.rtm,...          % whether to realign to mean 
            'interp', defs.estimate.interp,...    % interpolation type
            'wrap', defs.estimate.wrap,...        % wrap in (x) y (z)
            'weight', defs.estimate.weight,...    % weight some voxels more?
            'sep', defs.estimate.sep,...          % interpolation size (separation)
            'PW',''...                            % ??
            );
        
        % Flags to pass to routine to create resliced images
        % (spm_reslice)
        resFlags = struct(...
            'interp', defs.write.interp,...       % interpolation type
            'wrap', defs.write.wrap,...           % wrapping info (ignore...)
            'mask', defs.write.mask,...           % masking (see spm_reslice)
            'which', aap.tasklist.currenttask.settings.reslicewhich,...     % what images to reslice
            'mean', aap.tasklist.currenttask.settings.writemean);           % write mean image
        
        clear imgs;
        for j = 1:length(aap.acq_details.sessions) % 
            % get files in this directory
            P = aas_getimages(aap,i,j,aap.tasklist.currenttask.epiprefix);
            imgs(j) = {P};
        end
        
        % Run the realignment
        spm_realign(imgs, reaFlags);
        
        if (~isdeployed)
            % Save graphical output
            figure(spm_figure('FindWin'));
            print('-djpeg','-r50',fullfile(aas_getsubjpath(aap,i),'diagnostic_aamod_realign'));
        end;
        
        % Run the reslicing
        spm_reslice(imgs, resFlags);
        
case 'checkrequirements'
    
otherwise
    aas_log(aap,1,sprintf('Unknown task %s',task));
end;












