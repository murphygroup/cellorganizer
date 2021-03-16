function a = tz_objmajorangle(obj)
%TZ_OBJMAJORANGLE Major angle of an object.
%   A = TZ_OBJMAJORANGLE(OBJ) returns the major angle of OBJ, which is an
%   [object] or [point array].
%   
%   See also ML_MAJORANGLE

%   03-Feb-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

img = tz_obj2img(obj,[]);
a = ml_majorangle(img);