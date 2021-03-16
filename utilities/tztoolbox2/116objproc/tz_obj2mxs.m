function [medaxis,width] = tz_obj2mxs(obj)
%TZ_OBJ2MXS Convert an object into medial axis representation.
%   MEDAXIS = TZ_OBJ2MXS(OBJ) returns the medial axis of an [object] or 
%   [point array].
%   
%   [MEDAXIS,WIDTH] = TZ_OBJ2MXS(...) also returns the widths.
%   
%   See also

%   13-May-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

[medaxis,width] = tz_img2mxs(tz_obj2img(obj,[]));