function number_of_channels = get_number_of_channels( img_path )
% NUMBER_OF_CHANNELS = GET_NUMBER_OF_CHANNELS( IMG_PATH )
%
% Input
% * img (a valid OME.TIFF file)
%
% Output
% * number of channels

% Ivan E. Cao-Berg
%
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
    number_of_channels = [];
	getReport( err )
	return
end

number_of_channels = omeMeta.getChannelCount(0);
reader.close();
end%get_number_of_channels