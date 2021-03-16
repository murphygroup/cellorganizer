function nuclearModel2slml( nuclearShapeModel, nuclearTextureModel, fileID )
% NUCLEARMODEL2SLML Handles the writing of the nuclear model to an
% SLML file

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
% Created: May 8, 2007
% Last Update: August 27, 2010
%
% Copyright (C) 2010 Department of Computational Biology
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

fprintf( fileID, '\t\t\t\t\t\t%s\n', '<compartment name="nucleus">' );
fprintf( fileID, '\t\t\t\t\t\t\t%s\n', '<listOfPatterns>' );
fprintf( fileID, '\t\t\t\t\t\t\t\t%s\n', '<pattern>' );
fprintf( fileID, '\t\t\t\t\t\t\t\t\t%s\n', '<object>' );

%nuclear shape model
if isempty( nuclearShapeModel )
    warning( 'CellOrganizer: No nuclear shape model is present' );
	fprintf( fileID, '\t\t\t\t\t\t\t\t\t\t%s\n', '<shape>' );
	fprintf( fileID, '\t\t\t\t\t\t\t\t\t\t%s\n', '</shape>' );
else
    nuclearShapeModel2slml( nuclearShapeModel, fileID );
end

%nuclear texture model
%if isempty( nuclearTextureModel )
%    warning( 'CellOrganizer: No nuclear texture model is present' );
%	fprintf( fileID, '\t\t\t\t\t\t\t\t\t\t%s\n', '<texture>' );
%	fprintf( fileID, '\t\t\t\t\t\t\t\t\t\t%s\n', '</texture>' );
%else
%    nuclearTextureModel2slml( nuclearTextureModel, fileID );
%end

fprintf( fileID, '\t\t\t\t\t\t\t\t\t\t%s\n', '<position>' );
fprintf( fileID, '\t\t\t\t\t\t\t\t\t\t%s\n', '</position>' );
fprintf( fileID, '\t\t\t\t\t\t\t\t\t\t%s\n', '<frequency>' );
fprintf( fileID, '\t\t\t\t\t\t\t\t\t\t%s\n', '</frequency>' );

fprintf( fileID, '\t\t\t\t\t\t\t\t\t%s\n', '</object>' );
fprintf( fileID, '\t\t\t\t\t\t\t\t%s\n', '</pattern>' );
fprintf( fileID, '\t\t\t\t\t\t\t%s\n', '</listOfPatterns>' );
fprintf( fileID, '\t\t\t\t\t\t%s\n', '</compartment>' );