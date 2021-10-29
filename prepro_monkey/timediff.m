function [imdiff, g, slicediff] = timediff(imgs, flags)
% Analyses slice by slice variance across time series
% FORMAT [imdiff, g, slicediff] = timediff(imgs, flags)
%
% imgs   - string or cell or spm_vol list of images
% flags  - specify options; if contains:
%           m - create mean var image (vmean*), max slice var image
%               (vsmax*) and scan to scan variance image (vscmean*) 
%           v - create variance image for between each time point
%
% imdiff - mean variance between each image in time series
% g      - mean voxel signal intensity for each image
% slicediff - slice by slice variance between each image 
%
% Matthew Brett 17/7/00
  
if nargin < 1
  imgs = spm_get(Inf, '*.img', 'Select time series images');
end
if iscell(imgs)
  imgs = char(imgs);
end
if ischar(imgs)
  imgs = spm_vol(imgs);
end
if nargin < 2
  flags = 'm';
end

nimgs = size(imgs,1);
if isempty(nimgs) | nimgs < 2
  return
end
V1 = imgs(1);
Vr = imgs(2:end);

ndimgs = nimgs-1;
Hold = 0;

if any(flags == 'v') % create variance images
  for i = 1:ndimgs
    vVr(i) = makevol(Vr(i),'v',16); % float
  end
end
if any(flags == 'm') % mean /max variance 
  mVr = makevol(V1,'std_',16); 
  sVr = makevol(V1,'std_sc_',16);
  xVr = makevol(V1,'std_scmax_',16);
end

[xydim zno] = deal(V1.dim(1:2),V1.dim(3));

p1 = spm_read_vols(V1);
slicediff = zeros(ndimgs,zno);
g = zeros(ndimgs,1);
for z = 1:zno
  M = spm_matrix([0 0 z]);
  pr = p1(:,:,z);
  if any(flags == 'm')
    [mv sx2 sx mxvs]  = deal(zeros(size(pr)));
  end
  cmax = 0;
  for i = 1:ndimgs
    c = spm_slice_vol(Vr(i),M,xydim,Hold);
    v = (c - pr).^2;
    slicediff(i,z) = mean(v(:));
    g(i) = g(i) + mean(c(:));
    if slicediff(i,z)>cmax
      mxvs = v;
      cmax = slicediff(i,z);
    end
    pr = c;
    if any(flags == 'v')
      vVr(i) = spm_write_plane(vVr(i),v,z);
    end
    if any(flags == 'm')
      mv = mv + v;
      sx = sx + c;
      sx2 = sx2 + c.^2;
    end
  end
  if any(flags == 'm') % mean std etc
    sVr = spm_write_plane(sVr,sqrt(mv/(ndimgs-1)),z);
    xVr = spm_write_plane(xVr,sqrt(mxvs),z);
    mVr = spm_write_plane(mVr,sqrt((sx2-((sx.^2)/ndimgs))./(ndimgs-1)),z);    
  end      
end
% $$$ if ~findstr(spm('ver'), '2')
% $$$    spm_close_vol([vVr sVr xVr mVr]);
% $$$ end

g = [mean(p1(:)); g/zno];
imdiff = mean(slicediff')';

return

function Vo = makevol(Vi, prefix, datatype)
Vo = Vi;
fn = Vi.fname;
[p f e] = fileparts(fn);
Vo.fname = fullfile(p, [prefix f e]);
Vo.dim(4) = datatype;
% if findstr(spm('ver'), '5')
% 
% %if findstr(spm('ver'), '2')
   Vo = spm_create_vol_fix(Vo, 'noopen');
% else
%   Vo = spm_create_image(Vo);
% end
return
    



