function img2=tz_maskimg_3d(img,maskimg)
%TZ_MASKIMG_3D Crop 3D image by an mask.
%
% OVERVIEW:
%   Creates a maskimg with the same z dimension as img and uses it to mask
%   the original image.

% Tao Peng
%
% Copyright (C) 2005-2016 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
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

%D. Sullivan 6/6/13 added support for 3D masks and catch/override for
%                   multiple masks being passed at once.

nslice=size(img,3);
%D. Sullivan 6/6/13 allow for 3D masks with the same size as the original
%image
if size(maskimg,3)~=nslice
    %D. Sullivan 6/6/13 force the mask without the same number of
    %slices to be 2D to start. otherwise take the first non-zero slice
    useslice = find(squeeze(sum(sum(maskimg,1),2)),1);
    if size(maskimg,3)>1
        warning(['Number of mask slices > 1 and not equal to number of ',...
            'image slices Using first non-zero mask channel. ',...
            'If multiple are cells present, please save each as a ',...
            'separate file with its own mask.']);
    end
    
    img2 = zeros(size(img));
    %D. Sullivan 6/6/13 get rid of for loop in favor of repmat for speed.
    %make sure the mask is 0's and 1's (could have been 0,255 or something)
    maskimg = repmat(maskimg(:,:,useslice),[1,1,nslice])>0;
    img2 = img.*maskimg;
else
    %D. Sullivan 6/6/13 again, no need for a for loop
    %make sure the mask is 0's and 1's (could have been 0,255 or something)
    img2 = (maskimg>0).*img;
end

