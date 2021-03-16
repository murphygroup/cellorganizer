function x = tz_ppoints(x1,x2,n)
%TZ_PPOINTS Generates points with equal interval
%   X = TZ_PPOINTS(X1,X2,N) return a row vector of N points betwee X1 and
%   X2.
%   
%   See also

%   02-Nov-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

% Copyright (C) 2007  Murphy Lab
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

if nargin < 3
    error('Exactly 3 arguments are required')
end

if n==1
    x = x1;
    return;
end

if n==2
    x = [x1 x2];
    return;
end

if x1==x2
    x = zeros(1,n)+x1;
    return;
end

interval = (x2-x1)/(n-1);
x = x1:interval:x2;
x(end) = x2;