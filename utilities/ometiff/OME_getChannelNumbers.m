function channel_numbers = OME_getChannelNumbers(imagepaths, channel_name)
%OME_getChannelNumbers returns a numel(imagepaths)x1 cell array, with the OMEtif channel index corresponding to channel_name.
% If channel_name is numeric, every entry of the returned cell array is channel_name
%
%Inputs:
% imagepaths =  Nx1 cell array of strings
% channel_name =  string or integer corresponding to the desired channel in each OMEtif
%
%Outputs:
% Nx1 cell array of channel indices

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


channel_numbers = cell(numel(imagepaths),1);

if isnumeric(channel_name)
    channel_numbers = num2cell(repmat(num2str(channel_name),size(imagepaths,1),1));
    return
end
for cur_path_ind = 1:numel(imagepaths)
    reader = bfGetReader(imagepaths{cur_path_ind});
    omeMeta = reader.getMetadataStore();
    num_channels = omeMeta.getPixelsSizeC(0).getValue();
    channel_names = cell(num_channels,1);
    for chan = 1:num_channels
        channel_names{chan} = omeMeta.getChannelName(0,chan-1);
    end
    channel_numbers{cur_path_ind} = find([channel_names{:}] == channel_name);
    reader.close();
end
end