function [n,frac] = tz_borderpixelnum(obj,imgsize)
%TZ_BORDERPIXELNUM Count the pixels along image borders in an object.
%   N = TZ_BORDERPIXELNUM(OBJ,IMGSIZE) returns the number of pixels located
%   at the border of an image with [image size] IMGSIZE, for the [object]
%   or [point array] OBJ.
%   
%   [N,FRAC] = TZ_BORDERPIXELNUM(...) also returns the fraction of such
%   pixles.
%   
%   See also

%   29-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('Exactly 2 arguments are required')
end

n = sum( (obj(:,1)==1) | (obj(:,2)==1) | (obj(:,1)==imgsize(1)) | ...
    (obj(:,2)==imgsize(2)) );

frac = n/size(obj,1);