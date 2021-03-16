function y = tz_spgfilta(x,param)
%TZ_SPGFILTA Angular component of an oriented subband.
%   Y = TZ_SPGFILTA(X,PARAM) returns the value of the function G for each 
%   element in X, where H is at p.55 in the paper [1]. PARAM has two
%   fields:
%       'K' - total number of orientations
%       'k' - the kth orientation, must be between 0 and K-1
%
%   [1] J. Portilla and E.P. Simoncelli, A Parametric Texture Model Based 
%   on Joint Statistics of Complex Wavelet Coefficents. International
%   Journal of Computer Vision 40(1), 49-71, 2000
%
%   See also TZ_SPHIGHFILTR

%   13-Mar-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

seg1 = find(abs(x-pi*param.k/param.K)<pi/2);

alphak = 2^(param.K-1)*factorial(param.K-1)/ ...
    sqrt(param.K*factorial(2*(param.K-1)));
y = zeros(size(x));
y(seg1) = alphak*(cos(x(seg1)-pi*param.k/param.K)).^(param.K-1);