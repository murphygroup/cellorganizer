function [ cropimg, range ] = pc12_cropImg( img, pad)
%PC12_CROPIMG Automatically crops images
%
%Inputs:
% img = array of image data
% pad = amount of padding to use
%
%Outputs:
% cropimg = cropped image
% range = range of signal within the resulting 2D crop

%Author:Gregory Johnson Date:Unknown
% Copyright (C) 2007-2013  Murphy Lab
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

warning([mfilename ' is depricated. Please use ''cropImg.m.'''])

if ~exist('pad', 'var')
    pad = 5;
end

imsize = size(img);
mask = img > 0;

%Crop the image so we dont waste extra compute time
sumx = find(sum(sum(mask,3),1));
sumy = find(sum(sum(mask,3),2));

xstart = sumx(1)-pad;
if xstart < 1; xstart = 1; end

xend = sumx(end)+pad;
if xend > imsize(2); xend = imsize(2); end

ystart = sumy(1)-pad;
if ystart < 1; ystart = 1; end

yend = sumy(end)+pad; 
if yend > imsize(1); yend = imsize(1); end


%Crop to remove empty space
cropimg = img(ystart:yend, xstart:xend,:); 

range = [ystart yend xstart xend];

end

