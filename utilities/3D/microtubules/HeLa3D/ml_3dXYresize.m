function newim = ml_3dXYresize(oldim, xyr)
% Sept 8, 2009, Aabid Shariff
% This is modified from ml_3dimresize, to perform:
% 1. Bicubic interpolation instead of bilinear interpolation.
% 2. Not perform any operation on the Z

% NOTE: This modified version uses bicubic instead of bilinear. Hence, the image is very clean. Hence, this may not be good for active contour segmentation to give the boundary since a smoothed version will be better there.

% FUNCTION NEWIM = ML_3DXYRESIZE(OLDIM, XYR, ZR)
% Resize 3d images.  All 3 dimensions will be halved in 
% ml_3dimresize(oldim, 0.5, 0.5)
% Note: If the downsample ratio are all (1/int), use ml_downsize will be more
% efficient.  However, please note the difference of the expression.
% ml_3dimresize(im, xyr, zr) is equivalent to 
% ml_downsize(im, [1/xyr, 1/xyr, 1/zr])
% oldim : the image to be resized
% xyr : resize ratio on xy plane
% zr: resize ratio on z
% newim: the resize image

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

% Xiang Chen
% Apr. 19, 2004

for m = 1 : size(oldim, 3)
    newim(:,:, m) = imresize(oldim(:,:,m), xyr);
end
