function [img2,rect] = ml_imcrop(img,param)
%ML_IMCROP Crop an image.
%   IMG2 = ML_IMCROP(IMG,PARAM) returns an image that is cropped from the
%   image [IMG]. The way of cropping is specified by the structure PARAM.
%   PARAM has a field 'method', which could be one of the following values:
%       'rect' : rectangle
%           'topleft' - the coordinate of the top left corner. If topleft
%               is empty, it will be found automatically so that the center
%               of IMG2 is the center of the centroid of IMG.
%           'size' - size of the rectangle.
%     
%   [IMG2,RECT] = ML_IMCROP(...) also returns the crop region in RECT,
%   which is and 2x2 matrix with the topleft coordinate in the first row
%   and the bottomrigth coordinate in the second row.
%
%   See also

%   28-Jun-2006 Initial write T. Zhao
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

switch(param.method)
    case 'rect'
        if isempty(param.topleft)
            centroid = ml_imcentroid(double(img));
            param.topleft = round(centroid)-floor(param.size/2);
        end
        bottomright = param.topleft+param.size-1;
        img2 = img(param.topleft(1):bottomright(1), ...
            param.topleft(2):bottomright(2));
        rect = [param.topleft;bottomright];
    otherwise
        error('Unrecognized cropping method');
end
