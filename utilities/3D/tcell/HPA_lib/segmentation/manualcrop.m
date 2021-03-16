function regionimage = manualcrop(img)
% MASKS = MANUALCROP(IMG)
% Allows a user to perform manual segmentation
% Input: IMG is a 2D, uint8 image that
% Output: MASKS is a 2D, uint16 image showing the user defined regions
%
% Usage:
%   - The user clicks in the image to defne a region of interest. On
%   the last point, the user must double click to let the function know
%   he/she is done defining the region.
%   - The user is prompted to end ('yes') or continue ('enter').
%   The user can then keep defining additional regions of interest until
%   complete.
%
% Modified from define_regions.m by Xiang Chen
% Last edited by Justin Newberg, 6/22/2007


% Preprocessing image to enhance contrast
if isa(img,'uint8')
   img = img - min(img(:));
   img = img * (255/max(img(:)));
end

simg = size(img);
regionimage = uint16(zeros(simg(1),simg(2)));
MarkerNo = uint16(0);
R = 'no';

img(:,:,3) = 0;
dispimg(img,logical(regionimage));
while( ~strcmp(R,'yes'))
   MarkerNo = MarkerNo + 1;
   mask = roipoly;
   regionimage(mask) = MarkerNo;
   disp(['Region no.' sprintf('%i',MarkerNo) ' defined']);
   dispimg(img,logical(regionimage));
   R = input('Done yet (''yes'' to end or hit enter to continue)?','s');
end

% This helper function displays the image
function dispimg(I,ch2,cmap)
I(:,:,3) = 255*uint8(ch2>0);

imshow(I);
if exist('cmap','var')
   colormap(cmap);
end

return
