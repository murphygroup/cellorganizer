function slml = xml2slml( filename )
%XML2SLML Parses an XML into a Matlab structure the SLML instance

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
% Created: June 2, 2008
%
% Copyright (C) 2008 Center for Bioimage Informatics/Murphy Lab
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

if nargin ~= 1
    error('CellOrganizer: Wrong number of input arguments');
elseif ~isaFile( filename )
    error('CellOrganizer: The input argument is not an existing file');
elseif ~isSlmlInstanceValid( filename )
    error('CellOrganizer: Input argument is not a valid SLML instance');
else
    data = xml2struct( filename );
    slml.timestamp = datestr(now);
    if ~isempty(length(data.Attributes))
        for i=1:1:length(data.Attributes)
            slml = setfield(slml, data.Attributes(i).Name, ... 
                data.Attributes(i).Value );
        end
    end
end
end%slml2struct
