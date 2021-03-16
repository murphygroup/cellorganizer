function answer = show_RGB_image_galaxy_wrapper( file, red, green, blue )
%SHOW_RGB_IMAGE_GALAXY_WRAPPER Helper function that reads an image and saves
%an RGB png file composed of the selected channels.

% Ivan E. Cao-Berg (icaoberg@cmu.edu)
%
% Copyright (C) 2016 Murphy Lab
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

disp(['Attempting to load image ' file '.']);
img = tif2img( file );
if isempty( img )
    answer = false;
    warning('Unable to read image file.');
    return
end

number_of_channels = size(img,3);
disp( ['The image has ' num2str(number_of_channels) ' channels.' ] );

if red > number_of_channels
    answer = false;
    warning('RED channel index is greater than the number of available channels.');
    return
else
    disp( ['Setting RED channel to index ' num2str(red) '.'] );
end

if green > number_of_channels
    answer = false;
    warning('GREEN channel index is greater than the number of available channels.');
    return
else
    disp( ['Setting GREEN channel to index ' num2str(green) '.'] );
end

if blue > number_of_channels
    answer = false;
    warning('BLUE channel index is greater than the number of available channels.');
    return
else
    disp( ['Setting BLUE channel to index ' num2str(blue) '.'] );
end

disp('Attempting to save RGB image.')
try
    img = cat( 3, cat(3, img(:,:,red), img(:,:,green)),  img(:,:,blue));
    imwrite( img, 'output.png', 'png' )
catch the_error
    getReport( the_error )
    warning( 'Unable to save image to disk.');
    return
end

answer = true;
