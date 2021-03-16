function object2 = tz_rotateobj(object,theta)
%TZ_ROTATEOBJ Rotate object.
%   OBJECT2 = TZ_ROTATEOBJ(OBJECT,THETA) rotates an object by THETA
%   degrees. There is no artificial hole in the new object. OBJECT must
%   have 3 columns.
%   
%   See also tz_rotateobjs

%   13-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

objimg = ml_obj2img(object,[],{'2d','og'});
objimg2 = imrotate(objimg,theta);
[m,n,c] = find(objimg2);
object2 = [m,n,c];
