function nuclearShapeModel2slml( model, fileID )
% NUCLEARSHAPEMODEL2SLML Handles the writing of the nuclear model to a
% SLML instance

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
% Created: May 8, 2007
% Last Update: August 27, 2010
%
% Copyright (C) 2010  Department of Computational Biology
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

%shape model
fprintf( fileID, '\t\t\t\t\t\t\t\t\t\t%s\n', '<shape>' );
fprintf( fileID, '\t\t\t\t\t\t\t\t\t\t\t%s\n', '<listOfParameters>' );
write( fileID, 'parameter', model.name, 12, 'name', 'name', ...
    'constant', 'false', 'complex', 'false', 'type', 'string' );
	
%open medial axis
fprintf( fileID, '\t\t\t\t\t\t\t\t\t\t\t\t%s\n', [ '<parameter name="medaxis" ' ...
    'constant="false" complex="true">' ] );
writeParameter( fileID, 'name', model.medaxis.name, 13, 'false', 'false', ...
    'string' );
write( fileID, 'parameter', model.medaxis.constknot, 13, 'name', ...
    'constknot', 'constant', 'true', 'complex', 'false', 'type', ...
    'double' );
write( fileID, 'parameter', model.medaxis.nknots, 13, 'name', ...
    'nknots', 'constant', 'true', 'complex', 'false', 'type', 'double' );
	
fprintf( fileID, '\t\t\t\t\t\t\t\t\t\t\t\t\t%s\n', [ '<parameter name="stat" ' ...
    'constant="false" complex="true">' ] );
write( fileID, 'parameter', model.medaxis.stat.name, 14, 'name', ...
    'name', 'constant', 'false', 'complex', 'false', 'type', 'string' );
write( fileID, 'parameter', model.medaxis.stat.mu, 14, 'name', ...
    'mu', 'constant', 'true', 'complex', 'false', 'type', 'double' );
write( fileID, 'parameter', model.medaxis.stat.sigma, 14, 'name', ...
    'sigma', 'constant', 'true', 'complex', 'false', 'type', 'double' );
fprintf( fileID, '\t\t\t\t\t\t\t\t\t\t\t\t\t%s\n', '</parameter>' ); 	
	 
%close medial axis 
fprintf( fileID, '\t\t\t\t\t\t\t\t\t\t\t\t%s\n', '</parameter>' );

%width
fprintf( fileID, '\t\t\t\t\t\t\t\t\t\t\t\t%s\n%s\n', [ '<parameter name="width" ' ...
     'constant="false" complex="true">' ] );

% if strcmpi( model.width.stat.name, 'mvn' )%multivariate normal distribution
    % write( fileID, 'parameter', model.width.stat.name, 13, 'name', ...
        % 'name', 'constant', 'false', 'complex', 'false', 'type', 'string', ...
        % 'notes', 'Multivariate Normal Distribution' );
% else%generic model without notes
    % write( fileID, 'parameter', model.width.stat.name, 13, 'name', ...
        % 'name', 'constant', 'false', 'complex', 'false', 'type', 'string', ...
        % 'notes', 'Multivariate Normal Distribution' );
% end

% write( fileID, 'parameter', model.width.constknot, 9, 'name', ...
    % 'constknot', 'constant', 'true', 'complex', 'false', 'type', ...
    % 'double' );
% write( fileID, 'parameter', model.width.nknots, 9, 'name', ...
    % 'nknots', 'constant', 'true', 'complex', 'false', 'type', 'double', ...
    % 'notes', 'Number of knots');

% fprintf( fileID, '\t\t\t\t\t\t\t\t\t%s\n', [ '<parameter name="stat" ' ...
    % 'constant="false" complex="true">' ] );

% if strcmpi( model.width.stat.name, 'mvn' )%multivariate normal distribution
    % write( fileID, 'parameter', model.width.stat.name, 13, 'name', ...
        % 'name', 'constant', 'false', 'complex', 'false', 'type', 'string', ...
        % 'notes', 'Multivariate Normal Distribution' );
% else%generic model without notes
    % write( fileID, 'parameter', model.width.stat.name, 13, 'name', ...
        % 'name', 'constant', 'false', 'complex', 'false', 'type', 'string', ...
        % 'notes', 'Multivariate Normal Distribution' );
% end
    
 % write( fileID, 'parameter', model.width.stat.mu, 13, 'name', ...
     % 'mu', 'constant', 'true', 'complex', 'false', 'type', 'double', ...
     % 'notes', 'Location parameter' );
 
 % write( fileID, 'parameter', model.width.stat.sigma, 13, 'name', ...
     % 'sigma', 'constant', 'true', 'complex', 'false', 'type', 'double' );

%fprintf( fileID, '\t\t\t\t\t\t\t\t\t%s\n', '</parameter>' );    
fprintf( fileID, '\t\t\t\t\t\t\t\t\t\t\t\t%s\n', '</parameter>' );
fprintf( fileID, '\t\t\t\t\t\t\t\t\t\t\t%s\n', '</listOfParameters>' );
fprintf( fileID, '\t\t\t\t\t\t\t\t\t\t%s\n', '</shape>' );
end%nuclearShapeModel2slml