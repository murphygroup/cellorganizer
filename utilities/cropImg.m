function [ cropimg, croprange ] = cropImg( img, pad)
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


if ~exist('pad', 'var')
    pad = 5;
end

imsize = size(img);

imgdims = ndims(img);

croprange = [];
cropstr = [];

for i = 1:imgdims
    
    sumdims = 1:imgdims;
    sumdims(i) = [];
    
    img_tmp = img;
    
    for j = 1:length(sumdims)
        img_tmp = sum(img_tmp, sumdims(j));
    end
    
    img_tmp = squeeze(img_tmp);
    
    crop_start = find(img_tmp,1) - pad;
    if crop_start < 1
        crop_start = 1;
    end
    
    crop_stop = find(img_tmp,1,'last') + pad;
    if crop_stop > imsize(i)
        crop_stop = imsize(i);
    end
    
    croprange = [croprange crop_start crop_stop];
    
    cropstr = [cropstr num2str(crop_start) ':' num2str(crop_stop)];
    if i~=imgdims
        cropstr = [cropstr ','];
    end
end
    
cropimg = eval(['img(' cropstr ');']);

end
