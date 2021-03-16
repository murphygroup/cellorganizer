function tz_imshowmasked(img,maskpos)
%TZ_IMSHOWMASKED Show image after masking
%   TZ_IMSHOWMASKED(IMG,MASKPOS) shows the 2D image IMG, in which pixels
%   not at one of the positions in MASKPOS will be shown as black. MASKPOS
%   must have two column with the first column for X coordinates and the
%   other one for Y coordinates.

%   ??-???-???? Initial write T. Zhao
%   02-NOV-2004 Modified T. Zhao
%       - use sub2ind
%       - change function name tz_show_img --> tz_imshowmasked
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

imgsize=size(img);
img2=zeros(imgsize);
maskind=sub2ind(imgsize,maskpos(:,1),maskpos(:,2));
img2(maskind)=img(maskind);

figure
imshow(img2,[]);