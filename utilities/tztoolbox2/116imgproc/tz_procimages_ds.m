function tz_procimages_ds(imagedir,savedir,task)
%TZ_PROCIMAGES_DS Drug image batch processing.
%   TZ_PROCIMAGES_DS(IMAGEDIR,SAVEDIR,TASK)
%   
%   See also

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function tz_projimages_ds(imagedir,savedir,task)
%
%OVERVIEW:
%   batch processing of ds images
%PARAMETERS:
%   imagedir - rootdir
%   savedir - saving directory for output
%   task - projection,mask or feature calculation
%DESCRIPTION:
%   
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   02-NOV-2004 Modified TINGZ

load diroot.mat

switch task
case 'make_projs'
    tz_process3_ds([imagedir '/' DIDIR],[savedir '/' DPDIR],[savedir '/' DMDIR],task,'tif', 65535,1,'cell');
case 'make_masks'
    tz_process3_ds([imagedir '/' DPDIR],[savedir '/' DMDIR],[savedir '/' DMDIR],task,'tif', 65535,1,'cell');
case 'calc_feats'
    tz_calfeats_3d_dsm(imagedir,'cell',65535)
case 'new_projs'
    
end