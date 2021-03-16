function [pts,newsize] = tz_impolgrid(imgsize)
%TZ_IMPOLGRID Grid in polar coordinate system.
%   PTS = TZ_IMPOLGRID([M N]) returns a 2X(MXN) matrix. Each column of
%   the matrix is the coordinate in the polar coordinate system. The
%   first row is a vector of angles with unit radius and range [-pi,pi].
%   The second row is a vector of radius with the range [0,sqrt(M^2+N^2)].
%   
%   [PTS,NEWSIZE] = TZ_IMPOLGRID(...) also returns the size of the grid.
%
%   See also

%   10-Oct-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

thstep = 0.2;
rstep = 1;
maxlen = norm(imgsize/2);

[x,y] = meshgrid(-pi:thstep:pi,0:rstep:maxlen);
x = x';
y = y';
newsize = size(x);
pts = [x(:) y(:)]';