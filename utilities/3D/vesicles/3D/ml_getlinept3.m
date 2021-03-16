function pts=ml_getlinept3(s,a,len,z)
%ML_GETLINEPT3 Get coordinates of points on a 3D line segment.
%   ML_GETLINEPT3(S,A,LEN,Z) returns coordinates of points on a line segment
%   with starting point S (3d), angle A, length LEN, ending in slice Z.
%   

% March 29, 2012 R.F.Murphy created from ml_getlinept2 by T. Zhao
% April 2, 2012 I. Cao-Berg Fixed estimation of zvals using the appropiate 
%                           indices
%
% Copyright (C) 2012  Murphy Lab
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

if nargin < 4
    error('4 arguments are required');
end

a=mod(a,360);
len=round(len);

switch a
case 0
    t=[s(1)+len,s(2)];
case 90
    t=[s(1),s(2)+len];
case 180
    t=[s(1)-len,s(2)];
case 270
    t=[s(1),s(2)-len];
otherwise
    ra=a*pi/180;
    t=round([s(1)+cos(ra)*len,s(2)+sin(ra)*len]);
end

%icaoberg march 29, 2012
pts(:,1:2)=ml_getlinept(s(1:2),t);

%icaoberg march 30, 2012
zvals=double([0:size(pts,1)-1])*double(z-s(3))/size(pts,1);
pts(:,3) = round(zvals+s(3));
