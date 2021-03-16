function img2=tz_filtimg2featimg(img,alpha,u0)
%TZ_FILTIMG2FEATIMG Convert a filtered image to to a feature image.
%   IMG2 = TZ_FILTIMG2FEATIMG(IMG,ALPHA,U0)
%   
%   See also

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function img2=tz_filtimg2featimg(img)
%
%OVERVIEW:
%   convert a filtered image to to a feature image
%PARAMETERS:
%   img - input image
%   alpha - parameter of nonlinear transform
%   u0 - frequency of the gabor filter
%RETURN:
%   img2 - output image
%DESCRIPTION:
%   for gabor texture segmentation
%   unsupervised texture segmentation using gabor filers, A.K. Jain et.al.
%

if ~exist('u0','var')
    u0=1;
end

bsize=max(size(img));
psize=2^ceil(log2(bsize));

sigma=0.5*psize/u0;

img2=(1-exp(-2*alpha.*img))./(1+exp(-2*alpha.*img));

% wndsize=round(2*sigma);

% for i=-wndsize:wndsize
%     for j=-wndsize:wndsize
%         w(i+wndsize+1,j+wndsize+1)=exp(-(i^2+j^2)/sigma^2);
%     end
% end

% img2=conv2(img2,w/(2*wndsize+1)^2,'same');