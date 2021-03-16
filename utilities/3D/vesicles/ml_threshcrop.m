function mask = ml_threshcrop(image, cropimage, method, scale )
% ML_THRESHCROP generates a cropped, thresholded, cleaned  image
% ML_THRESHCROP(IMAGE, CROPIMAGE), where image is the IMAGE 
%    to be processed and CROPIMAGE is a binary mask defining
%    a region of interest.  Use CROPIMAGE=[] to process the 
%    entirety of IMAGE.  Thresholding is done BEFORE applying
%    the region of interest.  The image is cleaned using the 
%    majority operation of bwmorph.
%
% 06 Aug 98
%   - 10 Jan 1999 : renamed to reflect the order of operations
%   - 21 Apr 2014 : allowed input of scale so that users can multiply the
%   the threshold value by a factor

% Copyright (C) 2006-2014 Murphy Lab
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

% Written by Michael Boland
% $Id: ml_threshcrop.m,v 1.5 2006/06/27 13:33:47 tingz Exp $

if nargin<3
    method = 'nih';
    scale = 1;
end

if nargin<4
    scale = 1;
end

switch method
case 'nih'
    Iscaled = ml_nihscale(image) ;
    Timage = scale * ml_threshold(Iscaled) ;
%     Ithresh = im2bw(Iscaled, Timage) ;
    Ithresh = Iscaled > Timage;
case 'rc'
%     Iscaled = uint8(floor(ml_rcscale(image))) ;
%     Timage = scale * ml_rcthreshold(Iscaled);
% %     Ithresh = im2bw(Iscaled, Timage/255) ;
%     Ithresh = Iscaled > Timage/255;
%     
    Ithresh = image > ml_rcthreshold(image);
end

% mask = bwmorph(Ithresh, 'majority') ;
mask = Ithresh ;

%
% If the crop image exists, make all pixels outside the masked area
%   equal to 0.
%
if (~isempty(cropimage))
%     mask = roifilt2(0, mask, ~cropimage) ;
    mask = mask & cropimage;
end

function scaledimage = ml_nihscale(image)
% ml_nihscale(image)
% scales the pixel values of an image to make it like an nih image
% with 256 grey levels
%
% W. Dirks, 1998
%

% $Id: ml_threshcrop.m,v 1.5 2006/06/27 13:33:47 tingz Exp $

s = image * 253/(max(image(:))-min(image(:)))+1;
s = s/255;

scaledimage = s;

function scaledimage = ml_rcscale(image)
%function scaledimage = ml_rcscale(image)

s = (image-min(image(:))) * 255/(max(image(:))-min(image(:)));
scaledimage = s;
