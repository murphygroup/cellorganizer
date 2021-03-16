function pt=tz_rotate_2d(orgpt,a)
%TZ_ROTATE_2D 2D rotation.
%   PT = TZ_ROTATE_2D(ORGPT,A) returns a [point array] that is the rotation
%   [point array] ORGPT by A rad.

%   ??-???-???? Initial write T. Zhao
%   23-Mar-2005 Add comments T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

tm=[cos(a) sin(a) 0; -sin(a) cos(a) 0;0 0 1];

pt=tz_transform_2d(orgpt,tm);