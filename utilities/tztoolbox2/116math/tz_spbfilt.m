function y = tz_spbfilt(x,param)
%TZ_SPBFILT Oriented filter for steerable pyramid.
%   Y = TZ_SPBFILT(X,PARAM) returns the value of the function B for each 
%   row in X, where B is at p.55 in the paper [1]. X must be a nx2 matrix,
%   and the first column is r and the second column is theta. PARAM has two
%   fields:
%       'K' - total number of orientations
%       'k' - the kth orientation, must be between 0 and K-1
%
%   [1] J. Portilla and E.P. Simoncelli, A Parametric Texture Model Based 
%   on Joint Statistics of Complex Wavelet Coefficents. International
%   Journal of Computer Vision 40(1), 49-71, 2000
%
%   See also TZ_SPHIGHFILTR TZ_SPGFILTA

%   13-Mar-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

r = x(:,1);
theta = x(:,2);

y = tz_sphighfiltr(r).*tz_spgfilta(theta,param);