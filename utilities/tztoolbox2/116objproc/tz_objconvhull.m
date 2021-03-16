function hull = tz_objconvhull(obj)
%TZ_OBJCONVHULL
%   HULL = TZ_OBJCONVHULL(OBJ)
%   
%   See also

%   17-May-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

k = convhull(obj(:,2),obj(:,1));
[boxcorner,boxsize] = tz_boundbox(obj(:,1:2));
hull = roipoly(zeros(boxsize),obj(k,1),obj(k,2));
