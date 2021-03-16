function img = ml_mask( img, maskimg)

% CLEANIMG = ML_MASK( IMG, MASKIMG)
% 
% Masks IMG with MASKIMG. MASKIMG must be binary.  Support either 3D maskimg
% or 2D muskimg.  If Maskimg is 2D, it will be reproduced to the whole stack
% Xiang Chen, Aug 14, 2002

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


if (size(img, 3) ~= size(maskimg, 3)) & (size(maskimg, 3) ~= 1)
    error('MASKING must be either 2D or have same number of slices as IMG.');
end


% cleanimg = repmat( uint8(0), size(img));
% maskstack = repmat( uint8(0), size(img));
if (size(maskimg, 3) > 1)
    img(find(maskimg==0))=0;
else
    for m = 1 : size(img, 3)
        tmpimg=img(:,:,m);
        tmpimg(find(maskimg==0)) = 0;
        img(:,:,m) = tmpimg;
    end
end
% mask = find( maskstack);
% cleanimg( mask) = img( mask);
