function ds = ml_bwdist2(img1,img2)
%ML_BWDIST2 Find distances between two images.
%   DS = ML_BWDIST2(IMG1,IMG2) returns a nx3 matrix. The first 2 columns
%   are coordinates of white pixels in IMG2 and the last column contains
%   the correpsonding distance to the object in IMG1. IMG1 and IMG2 should
%   be [binary image]s.
%   
%   See also

%   16-Oct-2006 Initial write T. Zhao
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


if nargin < 2
    error('Exactly 2 arguments are required');
end

distimg = bwdist(img1);

[x,y] = find(img2>0);

ds = distimg(sub2ind(size(distimg),x,y));
ds = [x,y,ds];
