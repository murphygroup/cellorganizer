function [hs,maxdist] = ml_hightcode_3d(img,param)
%ML_HIGHTCODE_3D Estimate the height of the voxels in a 3D image.
%   HS = ML_HIGHTCODE_3D(IMG) returns the height based on the distance of a
%   refrence slice, which is the last slice. IMG should be a 3D [binary
%   image].
%   
%   HS = ML_HIGHTCODE_3D(IMG,PARAM) specifies how to calculate the heights
%   according to PARAM, which is a structure with the following fields:
%       'ref' - the index of reference slice
%       'res' - resolution
%   
%   [HS,MAXDIST] = ML_HIGHTCODE_3D(...) also returns the maximum distance.
%
%   See also

%   16-Oct-2006 Initial write T. Zhao
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


if nargin < 1
    error('1 or 2 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end


param = ml_initparam(param,struct('ref',size(img,3),'res',[]));

refimg = img(:,:,param.ref);

if all(refimg(:)==0)
    warning('The reference image is empty');
    hs = [];
    return;
end

distimg = bwdist(refimg);
maxdist = max(max(distimg(imfill(refimg,'hole')==1)));

hs = [];

for i=1:size(img,3)
    if sum(sum(img(:,:,i)))>0        
        [x,y] = find(img(:,:,i)>0);
        ds = distimg(sub2ind(size(distimg),x,y));
        %ds = [x,y,ds];
        ds = mean(ds);
        height = i-param.ref;
        hs = [hs;ds,zeros(size(ds,1),1)+height];
    end    
end

