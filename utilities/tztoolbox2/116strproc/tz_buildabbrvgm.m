function gm = tz_buildabbrvgm(s,param)
%TZ_BUILDABBRVGM
%   GM = TZ_BUILDABBRVGM(S)
%   
%   GM = TZ_BUILDABBRVGM(S,PARAM)
%   
%   See also

%   06-Feb-2007 Initial write T. Zhao
%   Copyright (c) 2007 Murphy Lab
%   Carnegie Mellon University
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation; either version 2 of the License,
%   or (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301, USA.
%   
%   For additional information visit http://murphylab.web.cmu.edu or
%   send email to murphy@cmu.edu


if nargin < 1
    error('1 or 2 arguments are required');
end

strlen = length(s);

gm = eye(strlen);
gm = [zeros(strlen,1) gm; zeros(1,strlen+1)];

for i=1:strlen-1
    for j=i+1:strlen
        [fullword,stat] = tz_expandabbrv(s(i:j));
        if stat==1
            gm(i,j+1) = 1;
        end
    end
end
