function answer = tif2ometiff( list_of_input_images, output_filename, parameters )
% TIF2OMETIFF Helper function that converts the output from output.tifimage into output.
%
% Make sure the output filename has the for <filename>.ome.tif
%
% Ulani Qi (uhq@andrew.cmu.edu), Ivan Cao-Berg (icaoberg@cmu.edu)
%
% Copyright (C) 2013-2016 Murphy Lab
% Computational Biology Department
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.o
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

% Ulani Qi (uhq@andrew.cmu.edu)
%
% Copyright (C) 2017 Murphy Lab
% Computational Biology Department
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.o
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

answer = false;
if isempty( list_of_input_images )
    warning( ['List of input files is empty. Exiting method.'  ] );
    return
end

% indicator variable for whether the number of rows and cols of the output
% ometiff has been standardized
found_rows_cols = 0;
for c=1:1:length(list_of_input_images)
    file = list_of_input_images{c};
    disp( ['Loading ' file] )
    try
        img = tif2img( file );
    catch err
        warning('Unable to load image. Exiting method.');
        getReport( err)
        return
    end
    
    disp( ['Image size is ' num2str(size(img))] );
    % this fix standardizes the number of rows and cols of the output
    % ometiff to follow that of the first tiff in the list_of_input_images
    if ~found_rows_cols
        [rows, cols, ~] = size(img);
        found_rows_cols = 1;
    end
    
    for z=1:1:size(img,3)
        %remember the indices mean xyzct
        omeimg(:,:,z,c,1) = imresize(img(:,:,z), [rows cols]);
    end
end

metadata = createMinimalOMEXMLMetadata(omeimg);

%metadata.PhysicalSizeX
if isfield( parameters, 'PhysicalSizeX' )
    pixelSize = ome.units.quantity.Length( ...
        java.lang.Double(parameters.PhysicalSizeX), ome.units.UNITS.MICROM);
    metadata.setPixelsPhysicalSizeX(pixelSize, 0);
else
    warning('PhysicalSizeX not set. Exiting method');
    return
end

%metadata.PhysicalSizeY
if isfield( parameters, 'PhysicalSizeY' )
    pixelSize = ome.units.quantity.Length( ...
        java.lang.Double(parameters.PhysicalSizeY), ome.units.UNITS.MICROM);
    metadata.setPixelsPhysicalSizeY(pixelSize, 0);
else
    warning('PhysicalSizeX not set. Exiting method');
    return
end

%metadata.PhysicalSizeZ
if isfield( parameters, 'PhysicalSizeZ' )
    pixelSize = ome.units.quantity.Length( ...
        java.lang.Double(parameters.PhysicalSizeZ), ome.units.UNITS.MICROM);
    metadata.setPixelsPhysicalSizeZ(pixelSize, 0);
else
    disp('PhysicalSizeZ not set. Assuming you are attempting to create a 2D image');
end

%metadata.list_of_channel_labels
if ~isfield( parameters, 'list_of_channel_labels')
    warning( 'Mandatory parameter list_of_channel_labels does not exist. Exiting method.' );
    return
elseif length(parameters.list_of_channel_labels) ~= length(list_of_input_images)
    warning('The number of input images must be equal to the number of labels. Exiting method.' );
    return
else
end

%metadata.list_of_channel_labels
for index=1:1:length(parameters.list_of_channel_labels)
    channel_name = parameters.list_of_channel_labels{index};
    channel_index = index-1;
    metadata.setChannelName( java.lang.String(channel_name), 0, channel_index )
end

bfsave( omeimg, output_filename, 'metadata', metadata, 'Compression', 'LZW' );
answer = true;