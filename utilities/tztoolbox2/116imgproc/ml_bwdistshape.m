function img2 = ml_bwdistshape(img,d)
%ML_BWDISTSHAPE Shrink or expand a shape in an image. 
%   IMG2 = ML_BWDISTSHAPE(IMG,D) returns a binary image containing a shape
%   that is zoomed from the shape in the image IMG, which is a binary image
%   and suppposed to contian perimeter pixels only. If D is greater than 0,
%   the shape will expaned to distance D. If D is less than 0, the shape will
%   shrink to distance -D.
%   
%   See also

%   17-Oct-2006 Initial write T. Zhao
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

if d==0
    img2 = img;
    return;
end

distimg = bwdist(img);
img = imfill(img,'hole');


if d>0
    distimg(img==1) = 0;
else
    distimg(img==0) = 0;
end
    
distimg(distimg<abs(d)) = 0;

img2 = bwperim(distimg);
