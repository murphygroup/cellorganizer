function list = get_list_of_function_handles_from_ometiff( image_path, ...
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

if ~exist( image_path )
    warning( ['Image file ' image_path ' does not exist.'] );
    return
end

numChannels = get_number_of_channels( image_path );
if channel_index < 1 || channel_index > numChannels
    error('Invalid channel index');
end

% Check channel index argument
number_of_rois = get_number_of_rois( image_path );

dimensionality = get_ometiff_dimensionality( image_path );

% if strcmp(dimensionality, '2D') == 1
if number_of_rois == 0 
    list{1} = @() get_roi_by_roi_index(image_path, channel_index, 0 );
end
for roiIndex=1:1:number_of_rois
    list{length(list)+1} = ...
        @() get_roi_by_roi_index(image_path, channel_index, roiIndex );
end
% else
%     
% end
end%get_list_of_function_handles_from_ometiff