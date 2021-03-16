function idx = tz_objinimg(obj,imgsize)
%TZ_OBJINIMG Obsolete. See ML_OBJINIMG.
%   IDX = TZ_OBJINIMG(OBJ,IMGSIZE) returns the indices of points in the
%   [object] or [point array] which are outside of an image with image size
%   IMGSIZE. 
%   
%   See also

%   15-May-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

error(tz_genmsg('of','tz_objinimg','ml_objinimg'));

if nargin < 2
    error('Exactly 2 arguments are required')
end

idx = find(obj(:,1)<1 | obj(:,2)<1 | ...
    obj(:,1)>imgsize(1) | obj(:,2)>imgsize(2));
