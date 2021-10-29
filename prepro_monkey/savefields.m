function savefields(fnam,p)

% Copied from the bottom of aamod_SegCoregNorm2d.m 2/21/17 TMD

if length(p)>1, error('Can''t save fields.'); end;
fn = fieldnames(p);
if numel(fn)==0, return; end;
for i=1:length(fn),
    eval([fn{i} '= p.' fn{i} ';']);
end;
if str2double(version('-release'))>=14,
    save(fnam,'-V6',fn{:});
else
    save(fnam,fn{:});
end;