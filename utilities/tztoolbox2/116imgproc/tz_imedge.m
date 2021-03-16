function img2=tz_imedge(img,cropimg,option)
%TZ_IMEDGE Obsolete. See ML_IMEDGE.
%   IMG2 = TZ_IMEDGE(IMG,CROPIMG,OPTION) returns an image that contains the
%   edges of objects in the image IMAGE, which will be masked first by
%   CROPIMG before edge detection. If CROPIMG is empty([]), there will be
%   no masking. OPTION specifies how to do edge detection:
%       'bw' - use BWPERIM directly.
%       'ce' - cell edge
%       'nu' - nuclear edge
%   This function is especially designed for extracting edges of cells and 
%   nucleus.
%
%   See also

%   27-SEP-2004 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

error(tz_genmsg('of','tz_imedge','ml_imedge'));

if nargin < 3
    error('Exactly 3 arguments are required')
end

img=double(img);

switch(option)
case 'bw'
    if ~isempty(cropimg)
        img(cropimg==0)=0;
    end
    img2=bwperim(img);
case 'ce'
    common = ml_imgcommonpixel(img);
    img2=bwmorph(img>common,'erode');
    if ~isempty(cropimg)
        img2(cropimg==0)=0;
    end
    img2 = bwperim(imfill(img2>0,'hole'));
%     medimg=medfilt2(img,[5 5]);
%     img2=bwperim(imfill(medimg>0,'hole'));
case 'nu'
    [img,maskimage]=ml_preprocess(double(img),cropimg,'ml','yesbgsub');
    img2=bwperim(imfill(img>0,'hole'));
end
    