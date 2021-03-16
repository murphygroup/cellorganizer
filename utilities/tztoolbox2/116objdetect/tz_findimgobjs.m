function objects=tz_findimgobjs(protimgfile,cropimgfile,way)

%function objects=tz_findimgobjs(protimgfile,cropimgfile,way)
%   
%OVERVIEW:
%   find objects in an image
%PARAMETERS:
%   protimgfile - image for object finding
%   cropimgfile - mask image
%   way - 'ml' or 'mb'
%RETURN:
%   objects - cell array of objects
%DESCRIPTION:
%
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   31-OCT-2004 Modified TINGZ
%       - add comments

im=mv_readimage(protimgfile);
cropimage=mv_readimage(cropimgfile);

[procimage,prot_maskimage]=ml_preprocess(im,cropimage,way,'yesbgsub');

objects=tz_findobjs(procimage);