function pts=ml_getline2pts(startpoint,endpoint,len)
%ML_GETLINE2PTS Get coordinates of points on a 3D line segment.
%   ML_GETLINEPT3(S,A,LEN,Z) returns coordinates of points on a line segment
%   with starting point START and ending point END.
%   

% April 3, 2012 R.F.Murphy created
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

if nargin < 2
    error('2 arguments are required');
end

diffs=double(endpoint)-double(startpoint);
[maxdif,idx]=max(abs(diffs));
incs=diffs/(maxdif+1);
len=uint16(max(len,maxdif+2));
offsets=repmat(double([0:1:len-1]),3,1).*repmat(incs',1,len);
pts=round(offsets+repmat(startpoint',1,len))';
