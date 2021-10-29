% This function extracts beta, st error and 90% CI
% for the voxel at the current location in the SPM
% results viewer (or at a location specified in
% the third argument).
% [beta,sterr,ci90] = extractSPMData(xSPM,SPM,[coords]);
% 11/6/2010 J Carlin
function [cbeta,SE,CI] = extractSPMData(SPM,coords, contrast)

hReg = findobj('Tag','hReg'); % get results figure handle

% Use current coordinates, if none were provided
if ~exist('coords','var')
	coords = spm_XYZreg('GetCoords',hReg);
end

% This is mostly pulled out of spm_graph - Original
% [xyz,i] = spm_XYZreg('NearestXYZ',coords,xSPM.XYZmm);
% spm_XYZreg('SetCoords',xyz,hReg);
% XYZ     = xSPM.XYZ(:,i); % coordinates

%-Parameter estimates:   beta = xX.pKX*xX.K*y;
%-Residual mean square: ResMS = sum(R.^2)/xX.trRV
%----------------------------------------------------------------------
%beta  = spm_get_data(SPM.Vbeta, XYZ);
%ResMS = spm_get_data(SPM.VResMS,XYZ);
%Bcov  = ResMS*SPM.xX.Bcov;

beta  = spm_get_data(SPM.Vbeta, coords);
ResMS = spm_get_data(SPM.VResMS,coords);
Bcov  = ResMS*SPM.xX.Bcov;

CI    = 1.6449;
% compute contrast of parameter estimates and 90% C.I.
%------------------------------------------------------------------
%Ic = xSPM.Ic; % Use current contrast
%contrast_names = {'simulation', 'perception', 'counting', 'simulation_over_counting', 'perception_over_counting', 'perception_simulation_over_counting', 'movement'};
contrast_names = struct2cell(SPM.xCon);
contrast_names = transpose(squeeze(contrast_names(1, 1, :)));
%contrast_names = {'simulation_button', 'perception_button', 'counting_button'};
Ic = find(strcmp(contrast_names, contrast));
cbeta = SPM.xCon(Ic).c'*beta;
SE    = sqrt(diag(SPM.xCon(Ic).c'*Bcov*SPM.xCon(Ic).c));
CI    = CI*SE;