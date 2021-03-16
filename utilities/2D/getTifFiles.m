function listOfFiles = getTifFiles( directory )
%GETTIFFILES Return a cell array containing the tif files in the directory.
%If the directory does not contain tif files it returns an empty structure.

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
% Created: September 17, 2007
% Last Update: March 2, 2008
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

if nargin == 0
    directory = pwd;
elseif nargin > 1
    error( 'CellOrganizer: Wrong number of input arguments' );
end

if ~isaDir( directory )
    listOfFiles = [];
    warning('CellOrganizer: Directory does not exist');
else
    listOfFiles = [];
    directoryStructure = dir( directory );
    for i=3:1:length(directoryStructure)
        information = imfinfo( [ directory filesep directoryStructure(i).name ] );
        if length(information) > 1
            imageFormat = information(1).Format;
        else
            imageFormat = information.Format;
        end

        if strcmpi( imageFormat, 'tif' )
            listOfFiles{length(listOfFiles)+1} = [ directory filesep directoryStructure(i).name ];
        end
    end
end%getTiffFiles