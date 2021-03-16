function pts = ml_imcoords(imgsize,scale,offset)
%ML_IMCOORDS Coordinates of all pixels in an image.
%   PTS = ML_IMCOORDS(IMGSIZE) returns a 3x(MxNxP) matrix if IMGSIZE is
%   [M,N,P]. The coordinates are obtained column by column.
%
%   PTS = ML_IMCOORDS(IMGSIZE,SCALE) will rescale the coordinates by
%   1/SCALE.
%
%   PTS = ML_IMCOORDS(IMGSIZE,SCALE,OFFSET) will move the coordinates by 
%   OFFSET after scaling.

%   11-Sep-2005 Initial write T. Zhao

% Copyright (C) 2007  Murphy Lab
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

% grj 6/30/2013 Modified for arbitrary dimensionality

if nargin < 1
    error('At least 1 argument is required')
end

ndims = length(imgsize);

if nargin < 2
    scale = 1;
end

if nargin < 3
    offset = zeros(1,ndims);
end

if length(scale) == 1
    scale = zeros(1,ndims) + scale;
end

coords = cell(ndims,1);
panels = cell(ndims,1);

for i = 1:length(imgsize)
    coords{i} = 1:imgsize(i);
end

[panels{:}] = ndgrid(coords{:});

for i = 1:length(panels)
    panels{i} = panels{i}(:)/scale(i) + offset(i);
end

pts = [panels{:}]';
