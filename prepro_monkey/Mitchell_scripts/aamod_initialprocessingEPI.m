% AA module  - Initial copying/unzipping/trimming/reeorenting of EPI 
% data from Oxford
% Daniel Mitchell 18/04/2012

function [aap,resp]=aamod_initialprocessingEPI(aap,task,i)

resp='';

switch task
    case 'domain'
        resp='subject';   % this module needs to be run once per subject
    case 'whentorun'
        resp='justonce';  % should this be run everytime or justonce?
    case 'description'
        resp='Run initalprocessing for EPIs';
    case 'summary'
        resp='Reformat raw FSL data';
    case 'report'
        
    case 'doit'
        subjpath=aas_getsubjpath(aap,i);
        
        % copy data to session and structural folders
        cwd=pwd;
        cd(subjpath);
        fprintf('\nCopying raw data...')
        estr='ep2dfmri.nii';
        epi=spm_select('list',subjpath,estr);
        
        sesspath=aas_getsesspath(aap,i,1);
        if ~exist(sesspath,'dir'), mkdir(sesspath); end
        epipath=fullfile(sesspath,epi);
        copyfile(epi,epipath);
        [junk Einfo]=aas_runfslcommand(aap,sprintf('fslinfo %s',epipath));
         
        epipath=regexprep(epipath,'.gz$','');
        epipathz=[epipath '.gz'];  
        
        if exist(epipathz,'file'),
            %fprintf('\nUnzipping EPI...')
            unix(sprintf('gunzip -f %s',epipathz));
        end   
        
        % check if EPI dimensions need reorienting (from sphinx position)
        S=textscan(Einfo,'%s\t%s');
        if str2double(S{2}{4})<str2double(S{2}{3})
            % slices already in 3rd dimension;
        elseif str2double(S{2}{4})>str2double(S{2}{3})
            % swap dimensions to match template (this sometimes leaves it unzipped)
            fprintf('\nReorienting EPI...')
            cmd=sprintf('fslswapdim %s x z y %s',epipath,epipath);
            aas_runfslcommand(aap,cmd);
            if exist(epipathz,'file'),
                %fprintf('\nUnzipping EPI...')
                unix(sprintf('gunzip -f %s',epipathz));
            end
        end  
        
        % split 4d
        fprintf('\nSplitting 4D epi file...')
        cd(fileparts(epipath))
        cmd=sprintf('fslsplit %s %s',epipath,[aap.acq_details.subjects(i).mriname '-0001-']);
        aas_runfslcommand(aap,cmd);
        movefile(epipath,fullfile(fileparts(epipath),'4Ddata'));
        % unzip everything
        %         fprintf('\nUnzipping files...')
        %         cmd=sprintf('gunzip -r %s',fileparts(strucpath));
        %         unix(cmd);
        cmd=sprintf('gunzip -rf %s',fileparts(epipath));
        unix(cmd);
        
        % reset orientations (with code from spm_image)
        fprintf('\nResetting orientations...')
        P=spm_select('FPlist',fileparts(epipath),'^mi.*.nii$');
        spm_progress_bar('Init',size(P,1),'Resetting orientations',...
            'Images Complete');
        for p=1:size(P,1),
            V    = spm_vol(deblank(P(p,:)));
            M    = V.mat;
            vox  = sqrt(sum(M(1:3,1:3).^2));
            if det(M(1:3,1:3))<0, vox(1) = -vox(1); end;
            orig = (V.dim(1:3)+1)/2;
            off  = -vox.*orig;
            M    = [vox(1) 0      0      off(1)
                0      vox(2) 0      off(2)
                0      0      vox(3) off(3)
                0      0      0      1];
            spm_get_space(P(p,:),M);
            spm_progress_bar('Set',i);
        end;
        spm_progress_bar('Clear');
        
        % Now move dummy scans to dummy_scans directory
        fprintf('\nMoving dummy scans...')
        dummypath=fullfile(sesspath,'dummy_scans');
        if (~length(dir(dummypath)))
            [s w]=aas_shell(['mkdir ' dummypath]);
            if (s)
                aas_log(aap,1,sprintf('Problem creating directory for dummy scans\n%s',dummypath));
            end;
        end;
        dummies=aas_getimages(aap,i,1,'',0,aap.acq_details.numdummies-1);
        for d=1:size(dummies,1)
            cmd=['mv ' dummies(d,:) ' ' dummypath];
            [s w]=aas_shell(cmd);
            if (s)
                aas_log(aap,1,sprintf('Problem moving dummy scan\n%s\nto\n%s\n',dummies(d,:),dummypath));
            end;
        end;
        
        cd(cwd);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;