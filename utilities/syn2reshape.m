function syn = syn2reshape( filename )
% SYN2RESHAPE Reshapes a 3D ometiff image into a 2D representation 
%
% List Of Input Arguments     Descriptions
% -----------------------     ------------
% img                         a 3D binary or realvalued image
%
% List Of Outputs     Descriptions
% ---------------     ------------
% img                 a 2D representation of the image 

% AUthor: Ivan E. Cao-Berg
% Copyright (C) 2013-2019 Murphy Lab
% Lane Center for Computational Biology
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

img = [];
try
    reader = bfGetReader( filename );
catch err
    disp('Failed to get image reader');
    getReport( err )
    return
end

try
    omeMeta = reader.getMetadataStore();
catch err
    disp('Failed to get metadata store');
    getReport( err )
    return
end

try
    size_x = omeMeta.getPixelsSizeX(0).getValue();
    size_y = omeMeta.getPixelsSizeY(0).getValue();
    size_z = omeMeta.getPixelsSizeZ(0).getValue();
    size_c = omeMeta.getPixelsSizeC(0).getValue();
catch err
    disp('Failed to extract image information');
    getReport( err )
    return
end

img = [];
for i=1:1:size_c
    temp = OME_loadchannel(filename, i);
    if size(temp,3) > 1
        temp = reshape( temp, size(temp, 1), [] );
    end
    img = [img; temp];
end

img = uint8(img);
reader.close();
end