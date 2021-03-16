function gf=tz_gauskernel(FWHM,Dim)
%TZ_GAUSKERNEL 2d Gaussian kernel - fairly continuous.
%   GF = TZ_GAUSKERNEL(FWHM,DIM)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

% 2d Gaussian kernel - fairly continuous

sig = FWHM/sqrt(8*log(2));

fineness = 0.1;
[x2d,y2d] = meshgrid(-(Dim(2)-1)/2:fineness:(Dim(2)-1)/2,...
		 -(Dim(1)-1)/2:fineness:(Dim(1)-1)/2);
gf    = exp(-(x2d.*x2d + y2d.*y2d)/(2*sig*sig));
gf    = gf/sum(sum(gf))/(fineness^2);
figure
colormap hsv
surf(x2d+Dim(1)/2,y2d+Dim(2)/2,gf);