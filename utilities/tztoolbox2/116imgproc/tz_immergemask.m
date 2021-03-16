function img2 = tz_immergemask(img,mask)
%TZ_IMMERGEMASK
%   IMG2 = TZ_IMMERGEMASK(IMG,MASK)
%   
%   See also

%   31-Oct-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

mask = bwperim(mask);
img2 = img;
img2(mask) = max(img(:));