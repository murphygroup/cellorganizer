function mixparas=tz_imggaussmix(img,mask,kmin,kmax)
%TZ_IMGGAUSSMIX Extimate guassian mixture for histogram of an image.
%   MIXPARAS = TZ_IMGGAUSSMIX(IMG,MASK,KMIN,KMAX)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function mixparas=tz_imggaussmix(img,mask)

if isempty(mask)
    data=img(img>0);
else
    data=img(mask==1);
end

[bestk,bestpp,bestmu,bestcov,dl,countf] = mixtures4(data(:)',kmin,kmax,0,1e-5,0)

mixparas=struct('k',bestk,'pp',bestpp,'mu',bestmu,'cov',bestcov);
