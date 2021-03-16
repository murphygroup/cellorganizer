function proj=tz_projimg_3d(img,mask,option,dim)
%TZ_PROJIMG_3D Project 3D matrix.
%   TZ_PROJIMG_3D(IMG,MASK,OPTION) returns a 2D matrix that is the 
%   projection of the 3D image IMG. MASK is a 2D inary image for masking.
%   OPTION specifies the way of projection:
%       'max' - pick the brightest pixel along projection axis
%       'mean' - take the mean of pixels along projection axis
%       'sum' - take sum of pixels along projection axis
%   Here the projection axis is Z axis. 
%   Notice: There is no big difference between option 'mean' and 'sum'.
%
%   TZ_PROJIMG_3D(IMG,MASK,OPTION,DIM) also specifies the projection
%   axis:
%       1 - X axis
%       2 - Y axis
%       3 - Z axis
%
%   See also TZ_PROJIMG_DIR TZ_PROCIMAGES_DS TZ_PROJIMG_MCF

%   17-Sep-2005 Initial write T. Zhao
%   ??-??-???? Initial write T. Zhao
%   18-AUG-2004 Modified T. Zhao
%       - add comments
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('3 or 4 arguments are required')
end

if ~exist('dim','var')
    dim = 3;
end

maskidx = find(mask==0);

if dim~=3
    for k=1:size(img,3)
        tmpimg = img(:,:,k);
        tmpimg(maskidx) = 0;
        img(:,:,k) = tmpimg;
    end
end

switch option
case 'max'
    proj=max(img,[],dim);
case 'mean'
    proj=mean(img,dim);
case 'sum'
    proj=sum(img,dim);
end

proj = squeeze(proj);

if dim==1
    proj = flipud(proj');
end

%For fast masking
if dim==3
    if ~isempty(mask)
        proj(maskidx)=0;
    end
end