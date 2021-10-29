
dbstop if error
warning off MATLAB:mir_warning_changing_catch

%% add paths
warning off all
addpath /imaging/local/software/AA/release-3.01/
aa_ver301
addpath(fileparts(fileparts(mfilename('fullpath'))),'-begin')
addpath(fileparts(mfilename('fullpath')),'-begin')
warning on all

fprintf('\nUsing %s.\n',spm('ver'))
if ~(strcmpi('SPM8',spm('ver')) || strcmpi('SPM12',spm('ver')))
    error('Please use SPM8')
end
spm('Defaults','FMRI');

%% create aap structure from parameters and tasklist files
aap=aarecipe('aap_parameters_defaults.xml','aap_tasklist_minimal.xml');
aap.schema.acq_details.subjects.tag.ATTRIBUTE=struct('desc','tag','ui','text');
aap.aap_beforeuserchanges.acq_details.subjects.tag=[];
aap.aap_beforeuserchanges.schema.acq_details.subjects.tag=[];
aap.aap_beforeuserchanges.schema.acq_details.subjects.tag.ATTRIBUTE=struct('desc','tag','ui','text');

%% Define study-specific parameters

% Directory for analysed data:
aap.acq_details.root='/imaging/dm01/AndrewBell/RestingState/spm8_FirstScans';
% note that since 11/03/14 we will not necessarily use the first scan per
% animal, but rather the last scan prior to any significant intervention
warning off all; 
addpath(aap.acq_details.root); 
warning on all

% Number of dummies:
aap.acq_details.numdummies=6; % used by Mars et al.

%% Add sessions
aap=aas_addsession(aap,'RestingStateBlock');

%% Add subjects, and EPI series numbers (ordered as aas_addsession above)
%% IF THIS CHANGES, CHANGE SUBSET FOR GROUP LEVEL BELOW
%%% Ignore for HPA project, but suitable controls otherwise:
aap=aas_addsubject(aap,'mi00198',{}); aap.acq_details.subjects(end).tag='Nev_4'; % Neville session 4
aap=aas_addsubject(aap,'mi00165',{}); aap.acq_details.subjects(end).tag='Nig_3'; % Nigel session 3 (1442 EPIs)
%%% Suitable controls for HPA project:
aap=aas_addsubject(aap,'mi00003',{}); aap.acq_details.subjects(end).tag='Pas_0'; % Passion (one bad slice)
aap=aas_addsubject(aap,'mi00019',{}); aap.acq_details.subjects(end).tag='Piz_0'; % Pizza
aap=aas_addsubject(aap,'mi00037',{}); aap.acq_details.subjects(end).tag='Pes_0'; % Pest 
aap=aas_addsubject(aap,'mi00039',{}); aap.acq_details.subjects(end).tag='Pap_0'; % Papa 
aap=aas_addsubject(aap,'mi00048',{}); aap.acq_details.subjects(end).tag='Oct_0'; % Octavia 
aap=aas_addsubject(aap,'mi00074',{}); aap.acq_details.subjects(end).tag='Pot_0'; % Pot
aap=aas_addsubject(aap,'mi00112',{}); aap.acq_details.subjects(end).tag='Spi_2'; % Spiro session 2 (ISO) 
aap=aas_addsubject(aap,'mi00114',{}); aap.acq_details.subjects(end).tag='Sau_2'; % Saul session 2 (ISO)
aap=aas_addsubject(aap,'mi00119',{}); aap.acq_details.subjects(end).tag='Sid_0'; % Sidd
aap=aas_addsubject(aap,'mi00123',{}); aap.acq_details.subjects(end).tag='Sil_2'; % Silo (odd, negative SeedToSeed matrix)
aap=aas_addsubject(aap,'mi00127',{}); aap.acq_details.subjects(end).tag='Rub_0'; % Rubin (872 EPIs)
aap=aas_addsubject(aap,'mi00137',{}); aap.acq_details.subjects(end).tag='Nip_0'; % Nippy (odd SeedToSeed matrix, with more connectivity on LHS)
aap=aas_addsubject(aap,'mi00157',{}); aap.acq_details.subjects(end).tag='Noo_0'; % Noodle
aap=aas_addsubject(aap,'mi00203',{}); aap.acq_details.subjects(end).tag='Seb_2'; % Seb check scan #? (942 EPIs)
aap=aas_addsubject(aap,'mi00214',{}); aap.acq_details.subjects(end).tag='Ran_2'; % Ranger check scan #?
aap=aas_addsubject(aap,'mi00222',{}); aap.acq_details.subjects(end).tag='Ste_2'; % Stevie check scan #? (No LR symmetry visible in SeedToSeed?)
aap=aas_addsubject(aap,'mi00224',{}); aap.acq_details.subjects(end).tag='Sco_2'; % Scottie check scan #?
aap=aas_addsubject(aap,'mi00253',{}); aap.acq_details.subjects(end).tag='Pop_0'; % Pop
aap=aas_addsubject(aap,'mi00255',{}); aap.acq_details.subjects(end).tag='Pud_0'; % Puddle
aap=aas_addsubject(aap,'mi00264',{}); aap.acq_details.subjects(end).tag='Rum_0'; % Rum
aap=aas_addsubject(aap,'mi00265',{}); aap.acq_details.subjects(end).tag='Rai_0'; % Raita (1000 EPIs)
% % from Jerome
aap=aas_addsubject(aap,'mi99076',{}); aap.acq_details.subjects(end).tag='Orv_0'; % Orvil 
aap=aas_addsubject(aap,'mi99080',{}); aap.acq_details.subjects(end).tag='Ors_0'; % Orson % 2 odd periods of reduced differential variance
aap=aas_addsubject(aap,'mi99083',{}); aap.acq_details.subjects(end).tag='Oli_0'; % Oliver % a slightly noisy slice
aap=aas_addsubject(aap,'mi99100',{}); aap.acq_details.subjects(end).tag='Ple_0'; % Pleco
aap=aas_addsubject(aap,'mi99105',{}); aap.acq_details.subjects(end).tag='Puc_0'; % Puck % odd regular little spikes
aap=aas_addsubject(aap,'mi99109',{}); aap.acq_details.subjects(end).tag='Puz_0'; % Puzzle 
aap=aas_addsubject(aap,'mi99111',{}); aap.acq_details.subjects(end).tag='Pri_0'; % Pringle 
aap=aas_addsubject(aap,'mi99115',{}); aap.acq_details.subjects(end).tag='Pug_0'; % Pugsley
aap=aas_addsubject(aap,'mi99162',{}); aap.acq_details.subjects(end).tag='Pil_0'; % Pilau % some spiking % check structural as it's already been trimmed so has unusual FOV
% (native space figures are chopped, but normalised brain looks fine)
aap=aas_addsubject(aap,'mi99166',{}); aap.acq_details.subjects(end).tag='Pep_0'; % Pepper 
 aap=aas_addsubject(aap,'mi99230',{}); aap.acq_details.subjects(end).tag='Orc_0'; % Orca
aap=aas_addsubject(aap,'mi99454',{}); aap.acq_details.subjects(end).tag='Ptt_0'; % Potto % 3 spikes % also has a 1.5 mm RS series

[aap.acq_details.subjects.structuralfn]=deal(aap.acq_details.subjects.mriname);

%% just record # WM and CSF components and iterations (hard-coded here for lack of a better place)
wmc=mean([1 2 1 1 1 2 1 1 1 1 1 1 1 1 1 6 1 1 1 1 2 1 4 1 1 1 1 1 1 1 1 1 1 1 1]);
csfc=mean([1 6 5 6 3 1 4 1 4 2 3 3 6 2 3 6 3 6 3 2 2 5 6 4 3 1 2 1 2 3 2 1 2 4 1]);
nit=mean([7 6 6 5 7 6 6 8 8 7 7 9 7 6 6 7 7 6 8 7 8 9 6 11 8 9 8 9 8 9 5 6 9 6 7]);

%% Some other settings
aap.directory_conventions.subject_directory_format=3; % don't change them
aap.directory_conventions.subject_filenames_format=0; % user specified
aap.directory_conventions.subject_filenames={aap.acq_details.subjects.mriname};
aap.directory_conventions.T1template='/imaging/dm01/AndrewBell/RhesusMacaqueAtlasTemplates/atlasfiles/112RM-SL_T1.nii';
aap.directory_conventions.rawdataafterconversionprefix='mi';
aap.options.autoidentifytmaps=false;
aap.options.autoidentifystructural=false;
aap.options.autoidentifystructural_average=false;
aap.options.autoidentifyfieldmaps=false;
aap.options.copystructuraltocentralstore=false;
aap.spm.defaults.normalise.write.vox=[2 2 2]; % resolution of normalised EPIs
aap.spm.defaults.normalise.write.bb=[-35 -30 -10; 35 55 45];
aap.timeouts.busy=60*24*3; % allow 3-day jobs in parallel mode

%% setup initial directory structure, then prompt for raw data if necessary
try
    aa_doprocessing(aap);
catch
    %empty={};
    needEPI={};
    needT1={};
    basedirs={'/imaging/dm01/AndrewBell/RestingState/spm_conn3',...
        '/imaging/ab03/HPA_Study',...
        '/imaging/ab03/_CBU_Scans',...
        '/imaging/ab03/Control_Dataset',...
        '/imaging/dm01/JeromeSallet/rs-fmri','/imaging/dm01/JeromeSallet/rs-fmri','/imaging/dm01/JeromeSallet/rs-fmri'};
    T1dir={'','structurals','structurals','structurals',...
        'structural','structural','structural'};
    T1filt={'^mi\d.*_structural','hr_\dMgFx.nii.gz$','hr_\dMgFx.nii.gz$','hr_\dMgFx.nii.gz$',...
        '^raw.nii','^structural.nii','mm(1001)?.nii'};
    EPIdir={'','nii_files','nii_files','nii_files',...
        'resting_state','resting_state','resting_state'};
    EPIfilt={'ep2d','ep2d.*fmri.nii','ep2d.*fmri.nii','ep2d.*fmri.nii',...
        'ep2d.*mri(_2mm)?.nii','ep2d.*mri(_2mm)?.nii','ep2d.*mri(_2mm)?.nii'};
    localdirs=cell(1,length(basedirs));
    for bd=1:length(basedirs)
        [junk, temp]=spm_select('list',basedirs{bd},'nofiles');
        localdirs{bd}=cellstr(temp);
    end
    for s=1:length(aap.acq_details.subjects)
        subdir=aas_getsubjpath(aap,s);
        temp=dir(subdir);
        [pth nam]=fileparts(subdir);
        % try to copy T1
        T1dest=spm_select('FPlist',subdir,[nam '_structural']);
        if exist(T1dest,'file')
            fprintf('\nFound %s',T1dest);
        else
            ok=false;
            for bd=1:length(basedirs)
                T1=spm_select('FPlist',fullfile(basedirs{bd},localdirs{bd}{strncmpi(nam,localdirs{bd},length(nam))},T1dir{bd}),T1filt{bd});
                %T1=spm_select('FPlist',fullfile(basedir,nam),['^' nam '_structural']);
                if size(T1,1)>1
                    fprintf('\nUsing 1st of %g:',size(T1,1))
                    T1=deblank(T1(1,:));
                end
                T1dest=fullfile(subdir,[nam '_structural' regexprep(T1,'[^\.]*(\..*)','$1')]);
                if exist(T1,'file') && ~exist(T1dest,'file')
                    copyfile(T1,T1dest);
                    fprintf('\nCopied %s',T1);
                    ok=true;
                    break
                end
            end
            if ~ok
                fprintf('\nDid not find %s',T1dest);
                needT1{end+1}=subdir;
            end
        end
        % check for EPI
        EPIdest=spm_select('FPlist',subdir,[nam '_ep2dfmri']);
         if exist(EPIdest,'file')
            fprintf('\nFound %s',EPIdest);
         else
             ok=false;
             for bd=1:length(basedirs)
                 EPI=spm_select('FPlist',fullfile(basedirs{bd},localdirs{bd}{strncmpi(nam,localdirs{bd},length(nam))},EPIdir{bd}),EPIfilt{bd});
                 if size(EPI,1)>1
                     clear sz
                     for w=1:size(EPI,1)
                         temp=dir(deblank(EPI(w,:)));
                        sz(w)=temp.bytes; 
                     end
                     [junk ind]=max(sz);
                     EPI=deblank(EPI(ind,:));
                 end
                 %T1=spm_select('FPlist',fullfile(basedir,nam),['^' nam '_structural']);
                 EPIdest=fullfile(subdir,[nam '_ep2dfmri' regexprep(EPI,'[^\.]*(\..*)','$1')]);
                 if exist(EPI,'file') && ~exist(EPIdest,'file')
                     copyfile(EPI,EPIdest);
                     fprintf('\nCopied %s',EPI);
                     ok=true;
                     break
                 end
             end
             if ~ok
                 fprintf('\nDid not find %s',EPIdest);
                 needEPI{end+1}=subdir;
             end
         end
    end
    if ~isempty(needEPI)
        fprintf('\nPlease copy EPI mi*_ep2dfmri.nii[.gz] into:\n')
        disp(char(needEPI))
        fprintf('\noptionally along with any physiological monitoring files, in format *_ISO.csv *_O2.csv *_CO2.csv\n')
    end
    if ~isempty(needT1)
        fprintf('\nPlease copy T1 mi*_structural.nii[.gz] into:\n')
        disp(char(needT1))
%         fprintf('\noptionally along with any T2, in format mi*_T2.nii[.gz]\n')
%         fprintf('\nand skull-stripped T1, in format mi*_ssT1.nii[.gz]\n')
    end
    if ~isempty(needEPI) || ~isempty(needT1)
        return
    end
end
clear temp empty needEPI needT1

%% Add tasks
aap=aas_addtask(aap,'aamod_initialprocessingEPI'); % Reorient, split 4D, unzip, and move dummy scans
aap=aas_addtask(aap,'aamod_tsdiffana','',{'[previous]'}); % Some time series diagnostics
aap=aas_addtask(aap,'aamod_realign','',{'aamod_initialprocessingEPI'});  % Realign EPIs to each other
aap=aas_addtask(aap,'aamod_tsdiffana','r',{'aamod_realign'}); % Repeat time series diagnostics
%%% there might be a problem with slice-time correction? so for now I'll skip it (see older version)
aap=aas_addtask(aap,'aamod_initialprocessingStruc'); % Reorient, trim, unzip
aap=aas_addtask(aap,'aamod_SegCoregNorm2d','',{'[previous]'}); % Coregister acquistions, tissue segmentation and normalisation to (McLaren) template
% note that above is not using realigned data (not obvious which is best in this case of minimal movement)
aap=aas_addtask(aap,'aamod_norm_write','',{'[previous]'});  % Write normalised (ie template space) EPIs
aap=aas_addtask(aap,'aamod_smooth','w',{'aamod_norm_write'}); % smooth images
aap.tasksettings.aamod_smooth(1).FWHM=3; % used by Mars et al. and Hutchinson et al.; Mantini used 2mm
aap.tasksettings.aamod_smooth(1).masks{1}='structurals/^wGM.nii';
aap.tasksettings.aamod_smooth(1).masks{2}='structurals/^wWM.nii';
aap.tasksettings.aamod_smooth(1).masks{3}='structurals/^wCSF.nii';
aap=aas_addtask(aap,'aamod_makeROIs','w',{'aamod_SegCoregNorm2d'}); % make/copy native and template space ROIs, and more conservative tissue masks
aap.tasksettings.aamod_makeROIs.templaterois='/imaging/dm01/AndrewBell/caret_brain/djm/Grown_wMacaque_CC11composite_andMD'; % individual ROIs
aap.tasksettings.aamod_makeROIs.warptonative=false;


analname='FP6MDiterateVEoverlappriors4_n35_SSSrerun'; %_n35
aap=aas_addtask(aap,'aamod_SPMconnIterateROIs2','s3w',{'aamod_makeROIs'}); 
aap.tasksettings.aamod_SPMconnIterateROIs2(end).analname=analname;
aap.tasksettings.aamod_SPMconnIterateROIs2(end).bpfilt=[0.0025 0.05];
aap.tasksettings.aamod_SPMconnIterateROIs2(end).windows=1;
aap.tasksettings.aamod_SPMconnIterateROIs2(end).confoundepiprefix='s3w'; % always use within-tissue-smoothed data for the confounds, to minimise modelling noise?
aap.tasksettings.aamod_SPMconnIterateROIs2(end).roipath=...
    {'/imaging/dm01/AndrewBell/caret_brain/djm/Grown_wMacaqueMD_LargestCluster/',...
    '/imaging/dm01/AndrewBell/caret_brain/djm/Grown_wMacaque.F99.LR.CC11composite_plusminus0p75/'};
% now loading unmasked ROIs and doing masking within module using file found with GMfilt
aap.tasksettings.aamod_SPMconnIterateROIs2(end).roifilt=...
  {'^MD_(MFGam|MFGp|FEF|preSMA|AI|IPS)_B.nii',...
     '_B'};
aap.tasksettings.aamod_SPMconnIterateROIs2(end).priors='_prior'; % set word, or suffix to roifilt, or single mask within which to use "growthzones"
aap.tasksettings.aamod_SPMconnIterateROIs2(end).tissuepath='structurals';
aap.tasksettings.aamod_SPMconnIterateROIs2(end).GMfilt='^t90rwc1mmi\d*_structural.nii$';
aap.tasksettings.aamod_SPMconnIterateROIs2(end).WMfilt='XrwWM.nii';
aap.tasksettings.aamod_SPMconnIterateROIs2(end).CSFfilt='XrwCSF.nii';
aap.tasksettings.aamod_SPMconnIterateROIs2(end).extraconfoundroifilt='DrawnSSS.*.nii';
aap.tasksettings.aamod_SPMconnIterateROIs2(end).extraconfoundroipath='/imaging/dm01/AndrewBell/RestingState';
aap.tasksettings.aamod_SPMconnIterateROIs2(end).componentspertissue=6.99; % this means the minimum of 6 or however many components explain 99% variance
aap.tasksettings.aamod_SPMconnIterateROIs2(end).tissueconfoundhistory=1;
aap.tasksettings.aamod_SPMconnIterateROIs2(end).motionconfoundhistory='satterthwaite';
aap.tasksettings.aamod_SPMconnIterateROIs2(end).globmove=1;
aap.tasksettings.aamod_SPMconnIterateROIs2(end).confoundsvdthresh=1; % i.e. none
aap.tasksettings.aamod_SPMconnIterateROIs2(end).doneflagsuffix=analname; %  

%% reset analysis steps?
% unix(sprintf('cd %s; rm -rvf */*/done*norm_write*',aap.acq_details.root))

%% Do processing
aa_doprocessing(aap);
%aa_doprocessing_parallel(aap);
%aa_doprocessing_parallel(aap,'continue');
%aa_doprocessing_cluster(aap,90,24); aa_doprocessing(aap); % gb, hrs; final doprocessing resaves aap with all subjects

%% split half
subs=aap.acq_details.subjects;

aap.acq_details.subjects=subs(1:2:end);
analname='FP6MDiterateVEoverlappriors4_n35_SSSrerun_odd'; 
aap=aas_addtask(aap,'aamod_SPMconnIterateROIs2','s3w',{'aamod_makeROIs'}); 
aap.tasksettings.aamod_SPMconnIterateROIs2(end)=aap.tasksettings.aamod_SPMconnIterateROIs2(end-1);
aap.tasksettings.aamod_SPMconnIterateROIs2(end).analname=analname;
aap.tasksettings.aamod_SPMconnIterateROIs2(end).doneflagsuffix=analname; %
aa_doprocessing(aap);

aap.acq_details.subjects=subs(2:2:end);
analname='FP6MDiterateVEoverlappriors4_n35_SSSrerun_even';
aap=aas_addtask(aap,'aamod_SPMconnIterateROIs2','s3w',{'aamod_makeROIs'}); 
aap.tasksettings.aamod_SPMconnIterateROIs2(end)=aap.tasksettings.aamod_SPMconnIterateROIs2(end-1);
aap.tasksettings.aamod_SPMconnIterateROIs2(end).analname=analname;
aap.tasksettings.aamod_SPMconnIterateROIs2(end).doneflagsuffix=analname; %
aa_doprocessing(aap);

%% display table of diagnostic images
addpath /imaging/dm01
addpath /imaging/dm01/MoreTools/
aap.models={'FP6MDiterateVEoverlappriors4_n35_SSS'};
aas_fmri_report3(aap)
