function img = generate_sphere( radius, options )
% GENERATE_SPHERE Simple function that draws a sphere inside a square
% matrix.
%
% Options
% radius

% Ivan E. Cao-Berg
%
% Copyright (C) 2016 Murphy Lab
% Computational Biology Department
% School of Computer Science
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
if nargin < 1 || nargin > 2
    warning('Wrong number of input arguments. Exiting method.' );
    return
end

if radius <= 0
    warning('Radius must me be a nonnegative number. Exiting method.');
    return
end

if nargin == 1
    options = [];
end

options = ml_initparam(options, ...
    struct('image_size', 512 ) );
options = ml_initparam(options, ...
    struct('x0', 0 ) );
options = ml_initparam(options, ...
    struct('y0', 0 ) );
options = ml_initparam(options, ...
    struct('z0', 0 ) );

%step1: generate sphere
n=linspace(-1,1,options.image_size);
[x,y,z]=ndgrid(n,n,n);
img=(sqrt((x-options.x0).^2 + (y-options.y0).^2 + (z-options.z0).^2)<=1);
img(find(img~=0))=255;
img=uint8(img);