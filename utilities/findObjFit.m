function [pos, rotangle] = findObjFit( regionObjs, regionBound, regionWeight, objErode, objbound, rotateAngle)
%Region - binary mask of area to place objects, to be eroded by objWeight
%regionWeight (optional) - pdf of object placement, independent of obj
%obj - 3d matrix of an object to be placed
%objWeight(optional) - the eroding element for region
%rotate angle

% grj 9/19/13 - improved object placement in 3D images

%
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


if ~exist('regionWeight', 'var') | isempty(regionWeight)
    regionWeight = regionObjs./sum(regionObjs(:));
end


if ~exist('objErode', 'var') | isempty(objErode)
    
    objErode = double(obj)./sum(double(obj(:)));
else
    objErode = double(objErode);
    objErode = objErode./sum(objErode(:));
end
if ~exist('rotateAngle', 'var') | isempty(rotateAngle)
    rotateAngle = 360;
end


nRotations = floor(360/rotateAngle);

inds = cell(1, nRotations);
angles = cell(1, nRotations);
weights = cell(1, nRotations);

%D. Sullivan - 5/28 making a copy of regionObjs seems unnecessary. Let's
%just use it 
% regionObjs = regionObjs;
regionsize = size(regionObjs);
if length(regionsize) == 3 && regionsize(3) >1
    regionObjs = padarray(regionObjs, [0,0,1]);
end
centroids = zeros(nRotations,3);
for i = 1:nRotations
    convObj = imrotate(objErode, rotateAngle*(i-1));
     %D. Sullivan 5/28/14 - adding padding to ensure no touching! 
%     convObj = imdilate(imrotate(objErode, rotateAngle*(i-1)),ones(2,2,2));
%     convObj = convObj./sum(convObj(:));
    %D. Sullivan and G. Johnson 5/22/14 - added corr to convnfft_fast to
    %use correlation instead of convolution so the objects aren't getting
    %flipped
    corr = 1;
    erodeRegion = convnfft_fast(regionObjs, convObj,corr);
    erodeRegion = erodeRegion.*regionObjs;
    if islogical(regionObjs)
        erodeRegion = erodeRegion > 0.999999;
    end
    
    if exist('regionBound', 'var') && exist('objbound', 'var')
        boundObj = imrotate(objbound, rotateAngle*(i-1));
        
        erodeBound = convnfft_fast(regionObjs, boundObj,corr);
        %D. Sullivan 5/28/14 - regionBound is a whole image only being used
        %as a logical flag. This seems very inefficient. 
        if islogical(regionBound)
            erodeBound = erodeBound > 0.999999;
        end
        
        erodeRegion = erodeRegion .* erodeBound;
    end
    
    %removes padding 
    if length(regionsize) == 3 & regionsize(3) >1
        erodeRegion = erodeRegion(:,:,2:end-1);
    end
    
    weight = regionWeight;
    weight(~erodeRegion) = 0;
    
    
    
    inds{i} = find(erodeRegion);
    angles{i} = ones(size(inds{i})) * i;
    weights{i} = weight(inds{i});



%     centroids(i,:) = round(getCentroid(convObj));
    centroid = round(size(convObj)/2);
    
    if length(centroid) == 2
        centroid = [centroid,1];
    end

    centroids(i,:) = centroid;

end
inds = vertcat(inds{:});
angles = vertcat(angles{:});
weights = vertcat(weights{:});

if isempty(inds)
    pos = [-1, -1, -1];
    return;
end

try
    ind = randsample(1:length(inds), 1, true, weights);
catch
    disp('asdf')
end

angle = angles(ind);

placementIndex = inds(ind);
%[y, x, z] = ind2sub(size(regionObjs), placementIndex);
[y, x, z] = ind2sub(size(erodeRegion), placementIndex);

pos = [y,x,z];

rotangle = rotateAngle*(angle-1);
% 
% objImg = placeObject(obj, targetSize, [y,x,z]);

end

function centroid = getCentroid(obj)
    objsize = size(obj);
    if length(objsize) == 2
        objsize = [objsize 1];
    end
    
    [x,y,z] = meshgrid(1:objsize(2), 1:objsize(1), 1:objsize(3));
    
    x1 = x.*obj;
    y1 = y.*obj;
    z1 = z.*obj;
    
    centroid = [sum(y1(:)) sum(x1(:)) sum(z1(:))]./sum(obj(:));
end