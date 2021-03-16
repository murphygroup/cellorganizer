function list = get_list_of_function_handles_from_wildcards( string, ...
    channel_index )
% GET_LIST_OF_FUNCTION_HANDLES_FROM_OMETIFF


% Ivan E. Cao-Berg (icaoberg@andrew.cmu.edu)
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

list = {};

if nargin ~= 2
    warning('Wrong number of input arguments');
    return
end

files = ml_ls(string);

if isempty( files )
    warning(['No files found using pattern ' string '.']);
    return
end

for index=1:1:length(files)
    image_path = files{index};
    number_of_rois = get_number_of_rois( image_path );
    number_of_channels = get_number_of_channels( image_path );
    if channel_index < 1 || channel_index > number_of_channels
        warning('Invalid channel index');
        return
    end
    for roiIndex=1:1:number_of_rois
        list{length(list)+1} = ...
            @() get_roi_by_roi_index(image_path, channel_index, roiIndex );
    end
end


end%get_list_of_function_handles_from_ometiff