function answer = show_projection_galaxy_wrapper( directory, channel )
%SHOW_PROJECTION_GALAXY_WRAPPER Helper function that reads an image and saves
%an RGB png file composed of the selected channel.

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

answer = false;
files = dir( [ directory filesep '*.tif' ] );
if isempty( files )
    warning(['No files found in directory ' directory '.' ]);
    return
end

if channel > length(files)
    warning( ['The channel selected (' num2str(channel) ')  is ' ...
        'greater than the number of available channels (' ...
        num2str(length(files)) '). Forcing channel index to 1.'] );
    channel = 1;
end

file = files( channel ).name;

try
    img = tif2img( file );
catch the_error
    warning('Unable to read image file.' );
    getReport( the_error );
    return;
end

img = img(:,:,channel );
options.method = 'sum';
img = im2projection( img, options );
img = uint8(img);

try
    imwrite( img, 'output.png' );
catch the_error
    warning( 'Unable to save image file to disk.' );
    getReport( the_error );
    return;
end

answer = true;