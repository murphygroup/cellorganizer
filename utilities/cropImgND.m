function [ cropimg, range ] = cropImgND( img, pad )
warning('This funcation is depricated and will be removed. Please use cropImg.m')

%CROPIMGND Automatically crops images
%
% Inputs:
% img = array of image data
% pad = amount of padding to use
%
% Outputs:
% cropimg = cropped image
% range = range of signal within the resulting 2D crop

% Author: Gregory Johnson
% Copyright (C) 2007-2016  Murphy Lab
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
    pad = 0;
end

imsize = size(img);
mask = img > 0;

ndims_img = ndims(img);
%Crop the image so we dont waste extra compute time

range = [];
cropstr = [];

for i = 1:ndims_img
    dims_to_sum = 1:ndims_img;
    dims_to_sum(i) = [];
    mask_tmp = mask;
    for j = 1:length(dims_to_sum)
        mask_tmp = sum(mask_tmp, dims_to_sum(j));
    end
    
    dim_proj_inds = find(mask_tmp);
    
    
    imgstart = dim_proj_inds(1)-pad;
    if imgstart < 1; imgstart = 1; end
    
    imgend = dim_proj_inds(end)+pad;
    if imgend > imsize(i); imgend = imsize(i); end
    
    
    range = [range, imgstart, imgend];
    cropstr = [cropstr, num2str( imgstart) ':'  num2str(imgend)];
    
    if i ~=ndims_img
        cropstr = [cropstr ','];
    end
end

cropimg = eval(['img(' cropstr ');']);

end

