function [segimage] = active3Dsegment(image,varargin)
%ACTIVE3DSEGMENT Active Contour 3D Segmentation
% First inputs must be either:
%   startSliceNo, endSliceNo - ints indicating start and end slices, or
%   startMask - a matrix the same size as image, 
% Then the following optional arugments
%   display - display 
%   alpha - the curvature alpha value
%   wholeImage - if true, the active contouring occurs on the whole image,
%   rather than per-slice

% Author: Aabid Shariff
%
% Copyright (C) 2012 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% March 15, 2012 R.F. Murphy Start with previous mask to speedup
% March 21, 2012 R.F. Murphy Add display param; remove disconnected pieces
% June 26, 2012 G.R. Johnson added optional start mask
% July 1, 2012 G.R. Johnson added option to perform segmentation on the
% whole image, rather than per-slice
% April 23, 2013 D. Sullivan Added iteration option to the 3D segmentation,
% default 1000
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

%Defaults
mask = [];

display  = 0;
alpha = 0.7;
sliceinc = 1;
wholeImage = 0;
%D. Sullivan 4/23/13 set maxIter default
maxIter = 1000;

%If the first option is a matrix
if numel(varargin{1}) > 1
    mask = varargin{1};
    
    inds = find(sum(sum(mask,1),2));
    startSlice = inds(1);
    stopSlice = size(mask,3);
    
    startMask = mask(:,:,startSlice);
    
    nextVar = 2;
    maskProvided = 1;
    
    mask = startMask;
else
    startSlice = varargin{1};
    stopSlice = varargin{2};
    nextVar = 3;
    maskProvided = 0;
    
    % decide which order to process slices in - better to start with biggest
    [X1,Y1]=find(image(:,:,startSlice));
    [X2,Y2]=find(image(:,:,stopSlice));
    areabox1 = (max(X1)-min(X1)+1)*(max(Y1)-min(Y1)+1);
    areabox2 = (max(X2)-min(X2)+1)*(max(Y2)-min(Y2)+1);
    sliceinc = 1; % assume start -> stop
    if areabox2 > areabox1
        tmp = startSlice;
        startSlice = stopSlice;
        stopSlice = tmp;
        sliceinc = -1; % reverse direction
    end
    
    startMask = [];
end

if length(varargin) >= nextVar && ~isempty(varargin{nextVar})
    display = varargin{nextVar};
end

nextVar = nextVar + 1;

if length(varargin) >= nextVar && ~isempty(varargin{nextVar})
    alpha = varargin{nextVar};
end

nextVar = nextVar + 1;

if length(varargin) >= nextVar && ~isempty(varargin{nextVar})
    wholeImage = varargin{nextVar};
end

%D. Sullivan 4/23/13 If maxIter is specified then set the value. otherwise
%use default of 1000
nextVar = nextVar + 1;

if length(varargin) >= nextVar && ~isempty(varargin{nextVar})
    maxIter = varargin{nextVar};
end


segimage = false(size(image));


if wholeImage
    mask = activcon_seg(image,mask,display, alpha);

    segimage=getBiggetsObj(mask);
else
    % use the mask from previous slice as the estimated mask for next slice
    for j = 1:2

        for I = startSlice:sliceinc:stopSlice
            %D. Sullivan 4/23/13 changed call to activcon_seg to include
            %maxIter parameter
%             mask = activcon_seg(image(:,:,I),mask,display, alpha);
            mask = activcon_seg(image(:,:,I),mask,display, alpha,maxIter);

            if I == startSlice && sum(segimage(:)) > 1
                startMask = mask;
            end

            segimage(:,:,I)=getBiggetsObj(mask);

        end

        startSlice = startSlice - sliceinc;

        mask = startMask;
        sliceinc = sliceinc * -1;

        if ~maskProvided
            break;
        else
            stopSlice = 1;
        end


    end
end
end

function mask = getBiggetsObj(mask)
    % pick the biggest object (clear the others)
    cc=bwconncomp(mask);
    if cc.NumObjects>1
        numPixels = cellfun(@numel,cc.PixelIdxList);
        [~,bigind] = max(numPixels);
        for i=1:cc.NumObjects
            if i~=bigind
                mask(cc.PixelIdxList{i})=0;
            end
        end
    end
    mask=imfill(mask, 'holes');
end
