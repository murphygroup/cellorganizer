function y=tz_rgb2gray(rgb,coef)
%TZ_RGB2GRAY Convert a rgb image to a gray-level image.
%   Y = TZ_RGB2GRAY(RGB,COEF)
%   
%   See also

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function y=tz_rgb2gray(rgb,coef)
%
%OVERVIEW:
%   convert a rgb image to a gray image
%PARAMETERS:
%   rgb - mxnx3 rgb image
%   coef - 1x3 converting coefficients
%RETURN:
%   y - mxn gray image
%DESCRIPTION:
%   y=c1*r+c2*b+c3*g
%
%HISTORY
%   24-OCT-2004 Initial write TINGZ

s=size(rgb);
r=rgb(:,:,1);
g=rgb(:,:,2);
b=rgb(:,:,3);

y=reshape([r(:),g(:),b(:)]*coef',s(1:2));
