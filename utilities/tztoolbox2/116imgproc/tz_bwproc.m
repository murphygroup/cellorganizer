function img2 = tz_bwproc(img,prep,kernel)
%TZ_BWPORC Process binary image.
%   IMG2 = TZ_BWPORC(IMG,PREP,KERNEL)
%   
%   See also

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function img2 = tz_bwproc(img,prep,kernel)
%OVERVIEW
%   process binary image
%PARAMETERS
%   img - input image
%   prep - preprocessing
%   kernel - convolution kernel
%RETURN
%   img2 - output image
%DESCRIPTION
%   
%HISTORY
%   08-May-2005 Initial write TINGZ
%SEE ALSO
%   

if ~isempty(prep)
    switch prep
    case 'pm4'
        img2=bwperim(img,4);
    case 'pm8'
        img2=bwperim(img,8);
    case 'dsm'
        img2=bwdist(img);
    case 'pms'
        img2=bwperim(img,8);
        img3=double(img)-double(img2);
    end
else
    img2=img;
end

if ~isempty(kernel)
    img2=conv2(img2,kernel,'same');
end

if strcmp(prep,'pms')
    subimg=conv2(img3,kernel,'same');
    img2=img2-subimg;
end