function y = tz_splowfiltr(x)
%TZ_SPLOWFILTR The radius of the low pass filter for steerable pyramid.
%   Y = TZ_SPLOWFILTR(X) returns the value of the function L for each 
%   elment in X, where L is at p.55 in the paper [1].
%
%   [1] J. Portilla and E.P. Simoncelli, A Parametric Texture Model Based 
%   on Joint Statistics of Complex Wavelet Coefficents. International
%   Journal of Computer Vision 40(1), 49-71, 2000
%   
%   See also

%   13-Mar-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

seg1 = find(x>pi/4 & x<pi/2);
seg2 = find(x<=pi/4);
seg3 = find(x>=pi/2);

y(seg1) = 2*cos(pi/2*log2(4*x(seg1)/pi));
y(seg2) = 2;
y(seg3) = 0;