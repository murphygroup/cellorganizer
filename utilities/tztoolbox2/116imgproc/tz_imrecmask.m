function mask = tz_imrecmask(img)
%TZ_IMRECMASK
%   MASK = TZ_IMRECMASK(IMG)
%   
%   See also

%   31-Oct-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

imshow(img);
pts = tz_fig2imcoord(round(ginput(2))');

mask = zeros(size(img));
mask(pts(1,1):pts(1,2),pts(2,1):pts(2,2))=1;
imshow(tz_immergemask(img,mask),[]);

