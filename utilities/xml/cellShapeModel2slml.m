function cellShapeModel2slml( model, fileID )
% CELLSHAPEMODEL2SLML Handles the writing of the cell membrane model to an
% SLML file

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
% Created: May 8, 2007
% Last Update: February 18, 2008
%
% Copyright (C) 2007  Murphy Lab
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
% For additional information visit http://murphylab.web.cmu.edu/software or
% send email to murphy@cmu.edu

fprintf( fileID, '\t\t\t\t%s\n', '<compartment name="cell membrane">' );

%object model
fprintf( fileID, '\t\t\t\t\t%s\n', '<object>' );
%shape model
fprintf( fileID, '\t\t\t\t\t\t%s\n', '<shape>' );
fprintf( fileID, '\t\t\t\t\t\t\t%s\n', '<listOfParameters>' );
write( fileID, 'parameter', model.name, 8, 'name', 'name', 'constant', ...
    'false', 'complex', 'false', 'type', 'string' );
write( fileID, 'parameter', model.startangle, 8, 'name', 'startangle', ...
    'constant', 'true', 'complex', 'false', 'type', 'string' );
write( fileID, 'parameter', model.anglestep, 8, 'name', 'anglestep', ...
    'constant', 'true', 'complex', 'false', 'type', 'string' );

fprintf( fileID, '\t\t\t\t\t\t\t\t%s\n', [ '<parameter name="stat" ' ...
    'constant="false" complex="true">' ] );
write( fileID, 'parameter', model.stat.name, 9, 'name', 'name', ...
    'constant', 'false', 'complex', 'false', 'type', 'string' );
write( fileID, 'parameter', model.stat.mu, 9, 'name', 'mu', ...
    'constant', 'true', 'complex', 'false', 'type', 'double' );
write( fileID, 'parameter', model.stat.sigma, 9, 'name', 'sigma', ...
    'constant', 'true', 'complex', 'false', 'type', 'double' );
fprintf( fileID, '\t\t\t\t\t\t\t\t\t%s\n', ['<parameter name=' ...
    '"transform" constant="false" complex="true">'] );
writeParameter( fileID, 'funname', model.stat.transform.funname, 10, ...
    'false', 'false', 'string', 'function' );
fprintf( fileID, '\t\t\t\t\t\t\t\t\t\t%s\n', ['<parameter name=' ...
    '"param" constant="false" complex="true">'] );
writeParameter( fileID, 'ncomp', model.stat.transform.param.ncomp, 11, ...
    'true', 'false', 'double' );
writeParameter( fileID, 'basevec', model.stat.transform.param.basevec, 11, ...
    'true', 'false', 'double', 'array' );
writeParameter( fileID, 'offset', model.stat.transform.param.offset, 11, ...
    'true', 'false', 'double', 'array' );

fprintf( fileID, '\t\t\t\t\t\t\t\t\t\t%s\n', '</parameter>' );
fprintf( fileID, '\t\t\t\t\t\t\t\t\t%s\n', '</parameter>' );
fprintf( fileID, '\t\t\t\t\t\t\t\t%s\n', '</parameter>' );

fprintf( fileID, '\t\t\t\t\t\t\t%s\n', '</listOfParameters>' );
fprintf( fileID, '\t\t\t\t\t\t%s\n', '</shape>' );
fprintf( fileID, '\t\t\t\t\t%s\n', '</object>' );
fprintf( fileID, '\t\t\t\t%s\n', '</compartment>' );
end%cellShapeModel2slml