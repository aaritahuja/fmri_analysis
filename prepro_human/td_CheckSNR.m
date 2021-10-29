function AllSNR = td_CheckSNR(MRI_path, subjnum, rundirs)

% Inputs:
% <subjectBoldPath> full path string to the location of the subject bold
%   directory
% <rundirs> cell array of the names of the run directories
% Output:
%  <AllSNR> mean SNR across 10 randomly selected images in each run
% Example:
% >> AllSNR = td_CheckSNR('/gpfs_home/dbasu1/data/pSEQ/subjects/102/bold', {'001' '002' '003' '004' '005' '006'})

% Original by Chris Gagne, updates by Brittany Ciullo
% Last modified 11/11/15 by Theresa M. Desrochers

subjectBoldPath = fullfile(MRI_path, num2str(subjnum), 'bold');

AllSNR = NaN(length(rundirs), 1);
for run = 1:length(rundirs)
  
    runPath = fullfile(subjectBoldPath, rundirs{run});
    
    % Find all the raw functional files
    % assumes they have 'f' as a prefix
    files = dir(fullfile(runPath,'f*.nii'));
    
    % Take 10 random indicies into the files
    randomidxs = randsample(length(files), 10);
    
    noise = NaN(length(randomidxs), 8);
    signal = NaN(length(randomidxs), 1);
    for ri = 1:length(randomidxs)
        
        %%%% ADD automask + mean/std
        % mean(x value you want, :, :) will give you all the y and z values for
        % this x slice
        
        % Construct the path to the random image
        filename = fullfile(runPath, files(randomidxs(ri)).name);
        
        % load image
        [epi,pixdim,rotate,dtype] = readnifti(filename);
        global_mean=mean(mean(mean(epi(:,:,:))));
        
        % calculate noise region (0:10,0:10,0:10)
        % TMD note - not sure why there are 8 calculations here...
%         noise(ri,1) = mean(mean(mean(epi(2:11,2:11,2:11))));
%         noise(ri,2) =mean(mean(mean(epi(end-11:end-1,end-11:end-1,end-11:end-1))));
%         noise(ri,3) =mean(mean(mean(epi(end-11:end-1,end-11:end-1,2:11))));
%         noise(ri,4) =mean(mean(mean(epi(end-11:end-1,2:11,end-11:end-1))));
%         noise(ri,5) =mean(mean(mean(epi(2:11,end-11:end-1,end-11:end-1))));
%         noise(ri,6) =mean(mean(mean(epi(2:11,2:11,end-11:end-1))));
%         noise(ri,7) =mean(mean(mean(epi(end-11:end-1,2:11,2:11))));
%         noise(ri,8) =mean(mean(mean(epi(2:11,end-11:end-1,2:11))));

noise(ri,1) = mean(mean(mean(epi(2:9,2:9,2:9))));
noise(ri,2) =mean(mean(mean(epi(end-9:end-1,end-9:end-1,end-9:end-1))));
noise(ri,3) =mean(mean(mean(epi(end-9:end-1,end-9:end-1,2:9))));
noise(ri,4) =mean(mean(mean(epi(end-9:end-1,2:9,end-9:end-1))));
noise(ri,5) =mean(mean(mean(epi(2:9,end-9:end-1,end-9:end-1))));
noise(ri,6) =mean(mean(mean(epi(2:9,2:9,end-9:end-1))));
noise(ri,7) =mean(mean(mean(epi(end-9:end-1,2:9,2:9))));
noise(ri,8) =mean(mean(mean(epi(2:9,end-9:end-1,2:9))));
        
        % noise(ran,rundir) = mean(mean(mean(epi(2:11,2:11,2:11))));
        
        % Calculate signal region
        % threshold brain using mean - 1 std. Rough estimate. Feel free to improve by using native mask etc.
        % Originally signal was kept across runs, but don't see a reason
        %   why - TMD
        %signal(ri,run) = mean(epi(epi>(.8*global_mean))); 
        signal(ri) = mean(epi(epi>(.8*global_mean)));
    end
    avg_snr = mean2(signal)/mean2(noise);
    
    AllSNR(run) = avg_snr;
    
end
