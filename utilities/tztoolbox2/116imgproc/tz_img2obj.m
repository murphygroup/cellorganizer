function obj = tz_img2obj(img,mask)
%TZ_IMG2OBJ Convert an image to an object.
%   OBJ = TZ_IMG2OBJ(IMG) returns an [object] from an [IMG] image. All
%   pixels with intensity above 0 in IMG will be considered as being in the
%   object.
%   
%   OBJ = TZ_IMG2OBJ(IMG,MASK) only takes pixels with value 1 in the
%   [binary image] MASK as an object. But if MASK is empty ([]), it is the
%   same as TZ_IMG2OBJ(IMG).
%   
%   See also

%   14-May-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('1 or 2 arguments are required')
end

if ~exist('mask','var')
    mask = [];
end

if isempty(mask)
    mask = img>0;
end

[r,c] = find(mask>0);
index = sub2ind(size(img),r,c);

obj = [r,c,img(index)];