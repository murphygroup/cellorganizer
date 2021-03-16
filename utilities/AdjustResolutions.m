function [img] = AdjustResolutions(img,initres,finalres,isbinary,ceilfloor)
%ADJUSTRESOLUTIONS resized a 3D image according to the input and output resolutions given
%
%Inputs:
% img = array of image data
% initres = initial resolution at which the img is being passed in
% finalres = final resolution at which the img is being returned

%Author: Devin Sullivan Summer 2013
% Copyright (C) 2013 Murphy Lab
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

if nargin<4 || isempty(isbinary)
    isbinary =1;
end

if nargin<5 || isempty(ceilfloor)
    ceilfloor = 'floor';
end

isblur = 1;

%since the cell is built conditionally on the nucleus, they should
%ALWAYS have the same resolution
% initres = param.resolution.cell;%e.g. [0.05,0.05,0.35]
% finalres = param.resolution.objects;
% maxval = max(img(:));
param.outputres=finalres;

%nucleus
switch lower(ceilfloor)
    case 'ceil'
        finalsize_x = ceil(initres(1)./finalres(1).*size(img,1));
        finalsize_y = ceil(initres(2)./finalres(2).*size(img,2));
        finalsize_z = ceil(initres(3)./finalres(3).*size(img,3));
    case 'floor'
        finalsize_x = floor(initres(1)./finalres(1).*size(img,1));
        finalsize_y = floor(initres(2)./finalres(2).*size(img,2));
        finalsize_z = floor(initres(3)./finalres(3).*size(img,3));
    case 'round'
        finalsize_x = round(initres(1)./finalres(1).*size(img,1));
        finalsize_y = round(initres(2)./finalres(2).*size(img,2));
        finalsize_z = round(initres(3)./finalres(3).*size(img,3));        
    otherwise
        error('Unrecognized rounding option for adjusting resolutions')
end

if isblur
    
end

img = imresize(img,[finalsize_x finalsize_y],'bilinear');
% img = imresize(img,[finalsize_x finalsize_y],'bilinear');
%need to resize the z
img = tp_stretch3d(img,finalsize_z);
% img = (img>0).*maxval;

if isbinary
    img = (img>0);
end

% [y,x,z] = ind2sub(size(img),find(img == 1));

% %cell
% %Note: the recalculation of final size should be unnecessary since the
% %cell and nucleus images should always be the same size, but since the
% %arithmatic is trivially fast I recalculate to deal with potential
% %weird situations in the future(e.g. if we need the nuclear image to be
% %a smaller object that we add in to the cell image for space)DPS
% finalsize_x = floor(initres(1)./finalres(1).*size(param.cell,1));
% finalsize_y = floor(initres(2)./finalres(2).*size(param.cell,2));
% finalsize_z = floor(initres(3)./finalres(3).*size(param.cell,3));
% param.cell = imresize(param.cell,[finalsize_x finalsize_y],'bilinear');
% %need to resize the z
% param.cell = tp_stretch3d(param.cell,finalsize_z);
% param.cell = param.cell>0;