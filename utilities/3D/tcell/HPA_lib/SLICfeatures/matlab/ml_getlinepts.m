function pts=ml_getlinepts(s,t)
%ML_GETLINPTS get all discrete points along a line
%   PTS=ML_GETLINEPTS(S,T) get points along a straight line staring from 
%   S to T. Both S and T should be a 1x2 vector.

% Copyright (C) 2006  Murphy Lab
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

d=t-s;

if all(d==0)
    pts=s;
    return;
end

if s(1)==t(1)
    pts=[zeros(abs(d(2))+1,1)+s(1) (s(2):sign(d(2)):t(2))'];
    return;
end

if s(2)==t(2)
    pts=[(s(1):sign(d(1)):t(1))' zeros(abs(d(1))+1,1)+s(2) ];
    return;
end

if abs(d(1))>abs(d(2))
    x=(s(1):sign(d(1)):t(1))';
    y=round(s(2)+(x-x(1))*d(2)/d(1));
else
    y=(s(2):sign(d(2)):t(2))';
    x=round(s(1)+(y-y(1))*d(1)/d(2));
end

pts=[x,y];