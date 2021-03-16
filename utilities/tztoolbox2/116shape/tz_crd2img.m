function img = tz_crd2img(pts,param)
%TZ_CRD2IMG Obsolete. See ML_CRD2IMG.
%   IMG = TZ_CRD2IMG(PTS) returns an image that contains the [curve] PTS.
%   
%   IMG = TZ_CRD2IMG(PTS,PARAM) specifies the parameters of convertion.
%   Currently it contains a field 'tz_obj2img' which has two subfields, 
%   'imgsize' and 'mode'. These two subfileds are parameters for 2nd and
%   3rd arguments for the function TZ_OBJ2IMG. 
%   
%   See also

%   31-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

error(tz_genmsg('of','tz_crd2img','ml_crd2img'));

if nargin < 1
    error('1 or 2 arguments are required')
end

if ~exist('param','var')
    param = struct([]);
end

obj2imgParameters.imgsize = [];
obj2imgParameters.mode = [];
param = ml_initparam(param,struct('tz_obj2img',obj2imgParameters));
pts = tz_showpts_2d(pts,'ln',0);
img = tz_obj2img(pts,param.tz_obj2img.imgsize,param.tz_obj2img.mode);

