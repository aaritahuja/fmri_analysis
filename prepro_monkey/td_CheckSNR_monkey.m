function AllSNR = td_CheckSNR_monkey(subjectBoldPath, rundirs)

% Inputs:
% <subjectBoldPath> full path string to the location of the subject bold
%   directory
% <rundirs> cell array of the names of the run directories
% Output:
%  <AllSNR> mean SNR across 10 randomly selected images in each run
% Example:
% >> AllSNR = td_CheckSNR('/mnt/sdd1/mSEQ/subjects/101/bold', {'001' '002' '003' '004' '005' '006'})

% Original by Chris Gagne, updates by Brittany Ciullo
% Last modified 5/5/15 by Theresa M. Desrochers

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
        % 8 points are the 8 corners of the cube
        noise(ri,1) = mean(mean(mean(epi(2:6,2:6,2:6))));
        noise(ri,2) = mean(mean(mean(epi(end-6:end-1,end-6:end-1,end-6:end-1))));
        noise(ri,3) = mean(mean(mean(epi(end-6:end-1,end-6:end-1,2:6))));
        noise(ri,4) = mean(mean(mean(epi(end-6:end-1,2:6,end-6:end-1))));
        noise(ri,5) = mean(mean(mean(epi(2:6,end-6:end-1,end-6:end-1))));
        noise(ri,6) = mean(mean(mean(epi(2:6,2:6,end-6:end-1))));
        noise(ri,7) = mean(mean(mean(epi(end-6:end-1,2:6,2:6))));
        noise(ri,8) = mean(mean(mean(epi(2:6,end-6:end-1,2:6))));
        
        
        % noise(ran,rundir) = mean(mean(mean(epi(2:5,2:5,2:7))));
        
        % Calculate signal region
        % threshold brain using mean - 1 std. Rough estimate. Feel free to improve by using native mask etc.
        signal(ri) = mean(epi(epi>(.8*global_mean))); 
    end
    avg_snr = roundn(mean2(signal)/mean2(noise), -2);
    
    AllSNR(run) = avg_snr;
    
end
