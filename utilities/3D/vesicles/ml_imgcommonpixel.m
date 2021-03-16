function commonpixel = ml_imgcommonpixel(image,threshold)
% ML_IMGCOMMONPIXEL finds the most common pixel value in the input image. 
% [COMMONPIXEL] = ML_IMGSCALE(IMAGE) 
% It only works properly for integer image.

% Author: Michael Boland
% February 23, 1999
%
% Copyright (C) 1999-2012 Murphy Lab
% Carnegie Mellon University
%
% 07 Jul 05 T. Zhao handle negative values
% 05 Nov 05 T. Zhao add threshold
% 15 Nov 11 G. Johnson changed calculation of min/max
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

%t+ Feb. 23, 2006
image = double(image);
%t++

%gj august 29, 2012
imagemin = min(image(:)) ;
imagemax = max(image(:)) ;

%Added by T. Zhao
if ~exist('threshold','var')
    threshold = imagemax;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if imagemin==imagemax
    commonpixel=imagemin;
    return;
end

if imagemin<0
    image=image-imagemin;
    imagemax = imagemax-imagemin ;
end

%
% Generate a histogram where the number of bins is equal to the 
%  maximum pixel value + 1 (0..maxpixel).
%

%Changed
% imghist = hist(image(:)/imagemax, imagemax+1);
%to
imghist = histc(image(:), floor(imagemin):ceil(threshold));
%GRJ 6/6/2012

% %added by T. Zhao
% imghist(ceil(threshold)+2:end) = 0;
% %%%%%%%%%%%%%%%%%%%

[hmax,ihmax] = max(imghist) ;

commonpixel = imagemin + ihmax - 1 ; 
% if imagemin<0
%     commonpixel=commonpixel+imagemin;
% end
