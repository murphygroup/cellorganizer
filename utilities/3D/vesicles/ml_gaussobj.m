function obj = ml_gaussobj( sigma, param )
%ML_GAUSSOBJ An object from Gaussian distribution.
%   OBJ = ML_GAUSSOBJ(SIGMA) returns an object that is extracted from a
%   2D Gaussian distribution which has covariance matrix SIGMA. The object
%   contains no less than 95% energy of the Gaussian distribution.
%   
%   See also

% 26-Jan-2006 T. Zhao
%
% Copyright (C) 2007-2012 Murphy Lab
% Carnegie Mellon University
%
% August 4, 2012 D. Sullivan added param structure to the method, param.imagesize
%       contains the size of the cell image which the objects will be added. It is
%       passed to the ml_gaussimg method
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

if nargin < 1
    error('Exactly 1 argument is required')
end

%devins 8/4/2012
img = ml_gaussimg( sigma, param );
if max(img(:))>1
    keyboard
end

y = ml_wquantile(img(:),0.9);
img(img<y) = 0;

imageSize = size(img);

idx = find(img>0);
[X,Y,Z] = ind2sub(imageSize,idx);
obj = [X,Y,Z,squeeze(img(idx))];
