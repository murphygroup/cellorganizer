function img2 = tz_imcart2pol(img)
%TZ_IMCART2POL Under construction
%   Transform an image from Cartesian to polar coordinates.
%   IMG2 = TZ_IMCART2POL(IMG) returns an image that is the transformation
%   of IMG in polar coordinate system.
%   
%   See also

%   10-Oct-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

cartpts = tz_imcoords(size(img));
[th,r] = cart2pol(cartpts(1,:),cartpts(2,:));

[polarpts,newsize] = tz_impolgrid(size(img));

zi = interp2(th,r,img(:)',polarpts(1,:),polarpts(2,:));

zi = reshape(zi,newsize(1),newsize(2));

