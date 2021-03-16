function roiid = get_roi_id( img_path, roi_index )
% ROIID = GET_ROI_ID( IMAGE_PATH, ROI_INDEX )
%
% Input
% * IMG_PATH ( path to a valid OME.TIFF file )
% * ROI_INDEX ( the 1-indexed index of a ROI in the OME.TIFF image )
%
% Output
% * The ROI ID string corresponding to the input ROI index in the input image. 

% Ulani Qi (uhq@andrew.cmu.edu)

% Copyright (C) 2018 Murphy Lab
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

try
	reader = bfGetReader(img_path);
	omeMeta = reader.getMetadataStore();
catch err
	warning('This method requires BioFormats for Matlab');
    roiid = '';
	getReport( err )
	return
end

% Input ROI index is 1-indexed, actual metadata is 0-indexed so minus 1
roiid = omeMeta.getROIID(roi_index-1);
reader.close();
end%get_roi_id