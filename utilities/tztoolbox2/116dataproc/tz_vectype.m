function flag = tz_vectype(x)
%TZ_VECTYPE Test if the given data is a numerical vector.
%   TZ_VECTYPE(X) returns 1 if X is a numerical row vector, 2 if X is a
%   numerical column vector and 0 if X is not a numerical vector. 
%   Notice: a scalar is not considered as a vector, so the returned value
%   is 0.
%   
%   See also

%   23-Aug-2006 Initial write T. Zhao
%   Copyright (c) 2006 Murphy Lab
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
    error('Exactly 1 argument is required');
end

if ~isnumeric(x)
    flag = 0;
    return;
end

if isempty(x)
    flag = 0;
    return;
end

if size(x,1)>1
    if size(x,2)>1
        flag = 0;
    else
        flag = 2;
    end
else
    if size(x,2)>1
        flag = 1;
    else
        flag = 0;
    end
end
