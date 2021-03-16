function img = tz_gaussimg(sigma)
%TZ_GAUSSIMG Obsolete. See ML_GAUSSIMG.
%   IMG = TZ_GAUSSIMG(SIGMA) returns an image that has intensities with
%   a 2D Gaussian distribution. SIGMA is the 2x2 covariance matrix of the 
%   Gaussian distribution.
%   
%   See also

%   26-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

error(tz_genmsg('of','tz_gaussimg','ml_gaussimg'));

if nargin < 1
    error('Exactly 1 arguments are required')
end

imageSize = round([6*sqrt(sigma(1,1)),6*sqrt(sigma(2,2))]);

x = tz_imcoords(imageSize,1,-round(imageSize/2))';
img = reshape(mvnpdf(x,[0 0],sigma),imageSize(1),imageSize(2));

