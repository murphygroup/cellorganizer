function img = ometiff2tile( filename )
% OMETIFF2RESHAPE Helper function that reshapes an OME.TIFF so we can
% display it on screen.
%
% This method will fail if the timepoint is greater than one.

% Author: Ivan E. Cao-Berg
%
% Copyright (C) 2017-2019 Murphy Lab
% Computational Biology Department
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

img = [];
if nargin == ~1
  warning('Wrong number of input arguments');
  return
end

reader = bfGetReader( filename );
omeMeta = reader.getMetadataStore();
size_x = omeMeta.getPixelsSizeX(0).getValue();
size_y = omeMeta.getPixelsSizeY(0).getValue();
size_z = omeMeta.getPixelsSizeZ(0).getValue();
img = bfopen( filename );
img = img{1}; %cell array
size_c = length(img) / size_z;

if size_z == 1
    img2 = [];
    for i=1:1:length(img)
        img2 = [img2; img{i}];
    end
    img = img2;
else
    img2 = [];
    for i=1:1:length(img)
        img2 = [img2, img{i}];
    end

    img=[ img2(:,1:(size(img2,2)/3)); ...
        img2(:,size(img2,2)/3+1:2*size(img2,2)/3); ...
        img2(:, 2*size(img2,2)/3+1:size(img2,2)) ];
end
end
