function [img] = OME_loadchannel(image_path,channel)
%OME_loadchannel returns the image corresponding to channel in an OMEtif. image_path can include the channel number, separated with a comma from actual path.
%
%Inputs:
% image_path = string
% channel = integer (optional)
%
%Outputs:
% img = MxNxZ array from image_path

%Author: Tim Majarian 7/28/16
% Copyright (C) 2016  Murphy Lab
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

if nargin == 1
    parts = strsplit(image_path, ' ,');
    image_path = parts{1};
    channel = str2double(strtrim(parts{2}));
elseif nargin > 2
    error('only 1 or 2 inputs work')
else
end

reader = bfGetReader(image_path);
omeMeta = reader.getMetadataStore();
% dimorder = omeMeta.getPixelsDimensionOrder(0).getValue();

try
    channel = str2double(strtrim(channel));
catch
end

if ischar(channel)
    channel = strtrim(channel);
    num_channels = omeMeta.getPixelsSizeC(0).getValue();
    channel_names = cell(num_channels,1);
    for chan =1:num_channels
        channel_names{chan} = omeMeta.getChannelName(0,chan-1);
    end
    try
        channel_num = find([channel_names{:}] == channel);
    catch
        disp('no channel names provided')
        img = [];
        return;
    end
else
    channel_num = channel;
end

x_size = omeMeta.getPixelsSizeX(0).getValue();
y_size = omeMeta.getPixelsSizeY(0).getValue();
z_size = omeMeta.getPixelsSizeZ(0).getValue();

img = zeros(y_size,x_size,z_size);
whole_img = tif2img(image_path);

for z = 1:z_size
    z_slice = reader.getIndex(z-1,channel_num-1,0)+1;
%     img(:,:,z) = bfGetPlane(reader, z_slice);
    img(:,:,z) = whole_img(:,:,z_slice);
end
reader.close();
end