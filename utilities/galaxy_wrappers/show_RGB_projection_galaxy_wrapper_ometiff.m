function answer = show_RGB_projection_galaxy_wrapper_ometiff( filename , red, green, blue )
%SHOW_RGB_PROJECTION_GALAXY_WRAPPER Helper function that reads an image and save
%an RGB png file composed of the selected channels.

% Author: Xin Lu
%
% Copyright (C) 2017 Murphy Lab
% Carnegie Mellon University
%
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
% For additional information visit http://murphylab.web.cmu.edu/ or
% send email to murphy@cmu.edu

answer = false;
%%
% img = uint8(img);
%%%

img = {};
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

if red > size_c
    warning( ['The selected red channel (' num2str(red) ')  is ' ...
        'greater than the number of available channels (' ...
        num2str(size_c) ').'] );
    return;
else
    disp( ['The selected red channel index is ' num2str(red) '.'] );
end
img{1} = OME_loadchannel(filename,red);

if green > size_c
    warning( ['The selected green channel (' num2str(green) ')  is ' ...
        'greater than the number of available channels (' ...
        num2str(size_c) ').'] );
    return;
else
    disp( ['The selected green channel index is ' num2str(green) '.'] );
end
img{2} = OME_loadchannel(filename,green);

if blue > size_c
    warning( ['The selected blue channel (' num2str(blue) ')  is ' ...
        'greater than the number of available channels (' ...
        num2str(size_c) ').'] );
    return;
else
    disp( ['The selected blue channel index is ' num2str(blue) '.'] );
end
img{3} = OME_loadchannel(filename,blue);

disp('Attempting to save RGB image.')
try
    img = im2projection_RGB( img );
    imwrite( img, 'output.png', 'png' )
catch the_error
    getReport( the_error )
    warning( 'Unable to save image to disk.');
    return
end

answer = true;
