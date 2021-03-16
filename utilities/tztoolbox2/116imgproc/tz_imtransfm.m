function [img2,offset,mask] = tz_imtransfm(img,h,method)
%TZ_IMTRANSFM Projective transformation of a 2D image.
%   IMG2 = TZ_IMTRANSFM(IMG,H) transform an image projectively by applying 
%   the transformation matrix H, which must be a 3x3 matrix. IMG could
%   be a 2D matrix for gray level image or a 3D matrix for color image.
%   
%   IMG2 = TZ_IMTRANSFM(IMG,H,METHOD) is not available now. It is expected
%   to support interpolation.
%   
%   [IMG2,OFFSET] = TZ_IMTRANSFM(...) also returns the offset of the 
%   orginal point ([0,0]) in the new image.
%
%   [IMG2,OFFSET,MASK] = TZ_IMTRANSFM(...) returns a mask for img2. Pixels
%   with value 0 in MASK means invalid points in IMG2.

%   11-Sep-2005 Initial write T. Zhao

if nargin < 2
    error('2 or 3 arguments are required')
end

if nargin < 3
    method = 'nearest';
end

if (size(img,1) <=1 | size(img,2) <= 2)
    error('Invalid image');
end

if any(size(h) ~= 3)
    error('H must be a 3x3 matrix');
end

if(det(h) == 0)
    error('The transformation matrix must not be singular');
end

invH = inv(h);
invH = invH/norm(invH);

originalCorners = [1 1 size(img,1) size(img,1); ...
    1 size(img,2) size(img,2) 1];

newCorners = tz_homoproj(originalCorners,h);
offset = min(newCorners,[],2)-1;
newSize = round(max(newCorners,[],2) - offset+1)';

offset = round(offset);

newCoords = tz_addrow( tz_imcoords(newSize)',offset' );
newCoords = newCoords';

mapCoords = tz_homoproj(newCoords,invH);

bgLevel = 0;%mean(img(:));
img2 = zeros(newSize)+bgLevel;

newCoords = tz_imcoords(newSize);
% mapCoords = round(mapCoords);
% mapCoords = tz_addcol(mapCoords,round(offset)-offset);
% offset = round(offset);
% offset = [0;0];

% invalidIndices = find(mapCoords(1,:)<=0 | mapCoords(1,:)>size(img,1) ...
%     | mapCoords(2,:)<=0 | mapCoords(2,:)>size(img,2));
% 
% mapCoords(:,invalidIndices) = [];
% newCoords(:,invalidIndices) = [];

if(size(img,3) == 1)
    [xgrid,ygrid] = meshgrid(1:size(img,1),1:size(img,2));
    img2(sub2ind(size(img2),newCoords(1,:),newCoords(2,:))) = ...
        interp2(xgrid,ygrid,double(img'),mapCoords(1,:),mapCoords(2,:),method);
%         img(sub2ind(size(img),mapCoords(1,:),mapCoords(2,:)));      
        
else
    img3 = img2;
    for i = 1:size(img,3)
        img4 = img3;
        singleChannel = img(:,:,i);
%         img4(sub2ind(newSize,newCoords(1,:),newCoords(2,:))) = ...
%             singleChannel(sub2ind(size(singleChannel), ...
%             mapCoords(1,:),mapCoords(2,:)));
        [xgrid,ygrid] = meshgrid(1:size(singleChannel,1), ...
            1:size(singleChannel,2));
        img4(sub2ind(size(img2),newCoords(1,:),newCoords(2,:))) = ...
            interp2(xgrid,ygrid,double(singleChannel'),mapCoords(1,:), ...
            mapCoords(2,:),method);

        img2(:,:,i) = img4;
    end
end

mask = ~isnan(img2);
img2(~mask) = 0;

