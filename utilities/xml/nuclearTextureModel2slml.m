function nuclearTextureModel2slml( model, fileID )
% NUCLEARTEXTUREMODEL2SLML Handles the writing of the nuclear texture model to an
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

if nargin ~= 2
    error( 'CellOrganizer: Wrong number of input arguments' );
elseif ~isa( model, 'struct' )
    error( 'CellOrganizer: Input parameter model must be a struct' );
elseif isempty( model )
    error( 'CellOrganizer: Input parameter model must not be empty' );
elseif ~isfield( model, 'name' )
    error( 'CellOrganizer: Cannot identify model because mandatory parameter name is not a field' );
else
    switch lower( model.name )
        case 'pyr'
            if ~isfield( model, 'pixelStats' ) || ...
                    ~isfield( model, 'pixelLPStats' ) || ...
                    ~isfield( model, 'autoCorrReal' ) || ...
                    ~isfield( model, 'autoCorrMag' ) || ...
                    ~isfield( model, 'magMeans' ) || ...
                    ~isfield( model, 'cousinMagCorr' ) || ...
                    ~isfield( model, 'parentMagCorr' ) || ...
                    ~isfield( model, 'cousinRealCorr' ) || ...
                    ~isfield( model, 'parentRealCorr' ) || ...
                    ~isfield( model, 'varianceHPR' )
                %texture model
                warning( 'CellOrganizer: One or several parameters of the nuclear texture model were missing. Writing an empty texture model.' );
                fprintf( fileID, '\t\t\t\t\t\t%s\n', '<texture>' );
                fprintf( fileID, '\t\t\t\t\t\t%s\n', '</texture>' );
            else
                %texture model
                fprintf( fileID, '\t\t\t\t\t\t%s\n', '<texture>' );
                fprintf( fileID, '\t\t\t\t\t\t\t%s\n', '<listOfParameters>' );

                %texture model name
                write( fileID, 'parameter', model.name, 8, 'name', 'name', ...
                    'constant', 'false', 'complex', 'false', 'type', 'string' );

                %pixelStats
                write( fileID, 'parameter', model.pixelStats, 8, 'name', 'pixelStats', ...
                    'constant', 'true', 'complex', 'false', 'type', 'double' );
                
                %pixelLPStats
                write( fileID, 'parameter', model.pixelLPStats, 8, 'name', 'pixelLPStats', ...
                    'constant', 'true', 'complex', 'false', 'type', 'double' );
                
                %autoCorrReal
                write( fileID, 'parameter', model.autoCorrReal, 8, 'name', 'autoCorrReal', ...
                    'constant', 'true', 'complex', 'false', 'type', 'double' );
                
                %magMeans
                write( fileID, 'parameter', model.magMeans, 8, 'name', 'magMeans', ...
                    'constant', 'true', 'complex', 'false', 'type', 'double' );
                
                %cousinMagCorr
                write( fileID, 'parameter', model.cousinMagCorr, 8, 'name', 'cousinMagCorr', ...
                    'constant', 'true', 'complex', 'false', 'type', 'double' );
                
                %parentMagCorr
                write( fileID, 'parameter', model.parentMagCorr, 8, 'name', 'parentMagCorr', ...
                    'constant', 'true', 'complex', 'false', 'type', 'double' );
                
                %cousinRealCorr
                write( fileID, 'parameter', model.cousinRealCorr, 8, 'name', 'cousinRealCorr', ...
                    'constant', 'true', 'complex', 'false', 'type', 'double' );

                %cousinRealCorr
                write( fileID, 'parameter', model.cousinRealCorr, 8, 'name', 'cousinRealCorr', ...
                    'constant', 'true', 'complex', 'false', 'type', 'double' );
                
                fprintf( fileID, '\t\t\t\t\t\t\t%s\n', '</listOfParameters>' );
                fprintf( fileID, '\t\t\t\t\t\t%s\n', '</texture>' );
            end
        otherwise
            error( 'CellOrganizer: Unknown or unsupported nuclear texture model. Please refer to documentation for more information.' );
    end
end
end%nuclearShapeModel2slml