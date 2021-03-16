function img = generate_ellipsoid( a, b, c, options )
%GENERATE_ELLIPSOID Simple helper that draw an ellipsoid in a 3D array.

% Ivan E. Cao-Berg
%
% Copyright (C) 2016 Murphy Lab
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

%step0: check input arguments
img = [];
if nargin < 3 || nargin > 4
    warning('Wrong number of input arguments. Exiting method.' );
    return
end

if a <= 0 || b <= 0 || c <= 0
    warning('Parameters a, b and c must me be a nonnegative number. Exiting method.');
    return
end

if nargin == 3
    options = [];
end

options = ml_initparam(options, ...
    struct('image_size', 512 ) );
options = ml_initparam(options, ...
    struct('number_of_points', options.image_size*2 ) );

options = ml_initparam(options, ...
    struct('x0', 0 ) );
options = ml_initparam(options, ...
    struct('y0', 0 ) );
options = ml_initparam(options, ...
    struct('z0', 0 ) );

%step1: make ellipsoid
n=linspace(-options.image_size,options.image_size,options.image_size);
[x,y,z]=ndgrid(n,n,n);
img=(sqrt((x-options.x0).^2/a^2 + (y-options.y0).^2/b^2 + (z-options.z0).^2/c^2)<=1);
img(find(img~=0)) = 255;
img = uint8(img);