function objects=tz_findmatobjs(protimgfile,way)

%function objects=tz_findmatobjs(protimgfile,way)
%   
%OVERVIEW:
%   find objects in a matlab data file
%PARAMETERS:
%   protimgfile - matlab data file containing the image
%   way - 'ml' or 'mb'
%RETURN:
%   objects - cell array
%DESCRIPTION:
%
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   31-OCT-2004 Modified TIINGZ
%       - add comments

load(protimgfile);

cropimage=ones(size(selimg));

[procimage,prot_maskimage]=ml_preprocess(double(selimg),cropimage,way,'yesbgsub');

objects=tz_findobjs(procimage);