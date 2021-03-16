function img2=ml_setimgptspixel(img,pts)

%ML_SETIMGPTSPIXEL sets values for specified pixels
%   IMG2=ML_SETIMGPTSPIXEL(IMG,PTS) change values of 
%   pixels in IMG. The positions of changed pixels
%   are defined by the first two columns of pts, for row
%   and colomn indices. The third column contain pixel values
%   of it exists. If PTS only has two columns, the pixels
%   are set to 1.

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

img2=img;
imgsize=size(img);
if size(pts,2)<3
    pts(:,3)=1;
end

pts(pts(:,1)<=0 | pts(:,1)>imgsize(1),:)=[];
pts(pts(:,2)<=0 | pts(:,2)>imgsize(2),:)=[];

img2(sub2ind(imgsize,pts(:,1),pts(:,2)))=pts(:,3);