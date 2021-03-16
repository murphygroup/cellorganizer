function [filename,selimg]=tz_findmajorslice(imgdir,ext,mask)
%TZ_FINDMAJORSLICE Find an image with most total intensity in a directory.
%   FILENAME = TZ_FINDMAJORSLICE(IMGDIR,EXT,MASK) returns the file name of
%   the major slice of the 3D image in the directory IMGDIR. Here the major
%   slice means the slice with most total intensity. All files with 
%   extension EXT under IMGDIR will be considered as slices of the 3D
%   image. The slices will be masked by the mask image MASK first before
%   intensity calculation.
%   
%   [FILENAME,SELIMG] = TZ_FINDMAJORSLICE(...) also returns the selected
%   slice, which is a masked image.
%   
%   See also

%   ??-???-???? Initial write T. Zhao
%   31-OCT-2004 Modified T. Zhao
%       - add comments
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('Exactly 3 arguments are required')
end

imgfiles=ml_dir([imgdir '/*.' ext]);
imgpath=[];
maxfluo=0;

for i=1:length(imgfiles)
    img=ml_readimage([imgdir '/' imgfiles{i}]);
    img(mask==0)=0;
    fluo=sum(img(:));
    if fluo>maxfluo
        selimg=img;
        maxfluo=fluo;
        filename=imgfiles{i};
    end
end
