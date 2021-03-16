function [center,mangle] = ml_objcenter(pts)
%ML_OBJCENTER Find the center of a contour.
%   CENTER = ML_OBJCENTER(EDGE) returns the center of a solid object 
%   specified by pts. If the number of
%   
%   [CENTER,MANGLE] = ML_OBJCENTER(...) also returns the major angle of
%   the object. The unit is radian (same as the function TZ_BWMAJORANGLE).
%   
%   See also ML_EDGECENTER

%   24-Apr-2005 Initial write  T. Zhao
%   27-Jan-2006 Modified T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

% Copyright (C) 2007  Murphy Lab
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

if nargin < 1
    error('Exactly 1 argument is required')
end

%need to find out when this if will be executed
if size(pts,2)==3
%   conpts=ml_showpts_2d(round(pts),'ln',1);
    edgeimg = ml_obj2img(pts,[]);
    objimg = imfill(edgeimg,'hole');
else
    %edgeimg = pts;
    objimg = pts;
end

s = regionprops(double(objimg),'Centroid');
center = s.Centroid;

% Find equator plane
%old form, slow changed by DPS 3/19/12
%for z = 1:size(objimg,3)
%    crossectionArea(z) = nnz(objimg(:,:,z));
%end
%[maxArea,equatorZ] = max(crossectionArea);
[maxArea,equatorZ] = max(sum(sum(objimg)));

% Find major angle
mangle = ml_bwmajorangle(objimg(:,:,equatorZ));

% center(1)=sum(pts(:,1).^2)/sum(pts(:,1));
% center(2)=sum(pts(:,2).^2)/sum(pts(:,2));
