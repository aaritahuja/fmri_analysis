% AA module - Initial copying/unzipping/trimming/reeorenting of structural 
% data from Oxford
% Daniel Mitchell 18/04/2012

function [aap,resp]=aamod_initialprocessingStruc(aap,task,i)

resp='';

switch task
    case 'domain'
        resp='subject';   % this module needs to be run once per subject
    case 'whentorun'
        resp='justonce';  % should this be run everytime or justonce?
    case 'description'
        resp='Run initalprocessing';
    case 'summary'
        resp='Reformat raw FSL data';
    case 'report'
        
    case 'doit'
        subjpath=aas_getsubjpath(aap,i);
        
        % copy data to session and structural folders
        cwd=pwd;
        cd(subjpath);
        fprintf('\nCopying raw data...')
        sstr='structural.nii';
        T2str='T2.nii';
        ssT1str='ssT1.nii';
        struc=spm_select('list',subjpath,sstr);
        T2=spm_select('list',subjpath,T2str);
        ssT1=spm_select('list',subjpath,ssT1str);
        
        strucpath=fullfile(subjpath,'structurals',struc);
        if ~exist(fileparts(strucpath),'dir'),
            mkdir(fileparts(strucpath));
        else
            delete(fullfile(fileparts(strucpath),'*'));
        end
        copyfile(struc,strucpath)
        [junk Sinfo]=aas_runfslcommand(aap,sprintf('fslinfo %s',strucpath));
        
        T2path=fullfile(subjpath,'structurals',T2);
        if ~isempty(T2)
            copyfile(T2,T2path)
            [junk T2info]=aas_runfslcommand(aap,sprintf('fslinfo %s',T2path));
        end
        
        ssT1path=fullfile(subjpath,'structurals',ssT1);
        if ~isempty(ssT1)
            copyfile(ssT1,ssT1path)
            [junk ssT1info]=aas_runfslcommand(aap,sprintf('fslinfo %s',ssT1path));
        end
        
        strucpath=regexprep(strucpath,'.gz$','');
        strucpathz=[strucpath '.gz'];
        T2path=regexprep(T2path,'.gz$','');
        T2pathz=[T2path '.gz'];
        ssT1path=regexprep(ssT1path,'.gz$','');
        ssT1pathz=[ssT1path '.gz'];
        
        if exist(strucpathz,'file'),
            %fprintf('\nUnzipping structural...')
            unix(sprintf('gunzip -f %s',strucpathz));
        end
        if exist(T2pathz,'file'),
            %fprintf('\nUnzipping T2...')
            unix(sprintf('gunzip -f %s',T2pathz));
        end
        if exist(ssT1pathz,'file'),
            %fprintf('\nUnzipping ssT1...')
            unix(sprintf('gunzip -f %s',ssT1pathz));
        end      
        
        % check if structural dimensions need reorienting (from sphinx position)
        S=textscan(Sinfo,'%s\t%s');
        if str2double(S{2}{4})<str2double(S{2}{3})
            % slices already in 3rd dimension; just unzip
        elseif str2double(S{2}{4})>str2double(S{2}{3})
            % swap dimensions to match template (this leaves it unzipped)
            fprintf('\nReorienting structural...')
            cmd=sprintf('fslswapdim %s x z y %s',strucpath,strucpath);
            aas_runfslcommand(aap,cmd);
            if exist(strucpathz,'file'),
                %fprintf('\nUnzipping structural...')
                unix(sprintf('gunzip -f %s',strucpathz));
            end
        end
        
        % trim structural if large FOV (this leaves it zipped)
        if str2double(S{2}{2})==512 % S{2}{2} should be large regardless of reorienting
            fprintf('\nTrimming structural...')
            cmd=sprintf('fslroi %s %s 128 256 128 256 0 -1',strucpath,strucpath);
            aas_runfslcommand(aap,cmd);
            if exist(strucpathz,'file'),
                %fprintf('\nUnzipping structural...')
                unix(sprintf('gunzip -f %s',strucpathz));
            end
        end
        
        % check if T2 dimensions need reorienting (from sphinx position)
        if ~isempty(T2)
            S=textscan(T2info,'%s\t%s');
            if str2double(S{2}{4})<str2double(S{2}{3})
                % slices already in 3rd dimension; just unzip
            elseif str2double(S{2}{4})>str2double(S{2}{3})
                % swap dimensions to match template (this leaves it unzipped)
                fprintf('\nReorienting T2...')
                cmd=sprintf('fslswapdim %s x z y %s',T2path,T2path);
                aas_runfslcommand(aap,cmd);
                if exist(T2pathz,'file'),
                    %fprintf('\nUnzipping T2...')
                    unix(sprintf('gunzip -f %s',T2pathz));
                end
            end
            
            % trim T2 if large FOV (this leaves it zipped)
            if str2double(S{2}{2})==512
                fprintf('\nTrimming T2...')
                cmd=sprintf('fslroi %s %s 128 256 128 256 0 -1',T2path,T2path);
                aas_runfslcommand(aap,cmd);
                if exist(T2pathz,'file'),
                    %fprintf('\nUnzipping T2...')
                    unix(sprintf('gunzip -f %s',T2pathz));
                end
            end
        end
        
        % check if ssT1 dimensions need reorienting (from sphinx position)
        if ~isempty(ssT1)
            S=textscan(ssT1info,'%s\t%s');
            if str2double(S{2}{4})<str2double(S{2}{3})
                % slices already in 3rd dimension; just unzip
            elseif str2double(S{2}{4})>str2double(S{2}{3})
                % swap dimensions to match template (this leaves it unzipped)
                fprintf('\nReorienting ssT1...')
                cmd=sprintf('fslswapdim %s x z y %s',ssT1path,ssT1path);
                aas_runfslcommand(aap,cmd);
                if exist(ssT1pathz,'file'),
                    %fprintf('\nUnzipping ssT1...')
                    unix(sprintf('gunzip -f %s',ssT1pathz));
                end
            end
            
            % trim ssT1 if large FOV (this leaves it zipped)
            if str2double(S{2}{2})==512
                fprintf('\nTrimming ssT1...')
                cmd=sprintf('fslroi %s %s 128 256 128 256 0 -1',ssT1path,ssT1path);
                aas_runfslcommand(aap,cmd);
                if exist(ssT1pathz,'file'),
                    %fprintf('\nUnzipping ssT1...')
                    unix(sprintf('gunzip -f %s',ssT1pathz));
                end
            end
        end     
        
        % reset orientations (with code from spm_image)
        fprintf('\nResetting orientations...')
        P=spm_select('FPlist',fileparts(strucpath),'^.*.nii$');
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
        
        cd(cwd);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;