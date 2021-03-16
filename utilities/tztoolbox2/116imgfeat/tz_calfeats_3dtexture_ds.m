function tz_calfeats_3dtexture_ds(rootdir,level,threshold)
%TZ_CALFEATS_3DTEXTURE_DS Obsolete.

%function tz_calfeats_3dtexture_ds(rootdir,level,threshold)
%
%OVERVIEW:
%   Calculate 3d texture features for DS tree
%PARAMETERS:
%   rootdir - root directory
%   level -  cell-level or object-level
%   threshold - loading threshold
%DESCRIPTION:
%   the image directory is given by tz_constant   
% 
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   26-MAR-2003 Modified TINGZ

error(tz_genmsg('of','tz_calfeats_3dtexture_ds'));

load diroot.mat
DFDIR='features_texture';
task='calc_feats_texture'

[features,names] = tz_process_all_3([rootdir '/' DIDIR],[rootdir '/' DMDIR],[rootdir '/' DFDIR],task,'tif',threshold,1,level);
save(['drug_' level '_' task '_feats_3.mat'],'features','names');

