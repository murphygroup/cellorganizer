function model = proteinModel2slml( model, compartmentName, ...
    proteinName, fileID )
% PROTEINMODEL2SLML Helper function that parses the protein model to SLML

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
% Last Update: February 15, 2008
%
% Copyright (C) 2008  Murphy Lab
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

proteinClass = getClass( compartmentName );
if isempty( compartmentName )
    compartment = '<compartment ';
else
    compartment = ['<comparment name="' compartmentName '" ' ];
end

if isempty( proteinName )
    protein = '';
else
    protein = ['protein="' proteinName '" '];
end

fprintf( fileID, '\t\t\t\t%s\n', [ compartment protein 'class="' proteinClass '">' ] );

writeParameter( fileID, 'name', model.name, 5, 'false', 'false', 'string' );

%object model
fprintf( fileID, '\t\t\t\t\t%s\n', '<object>' );
writeParameter( fileID, 'name', model.objectModel.name, 6, 'false', 'false', 'string' );
writeParameter( fileID, 'covartype', 'spherical', 6, 'false', 'false', ...
    'string' ); 
fprintf( fileID, '\t\t\t\t\t\t%s\n', [ '<parameter name="stat" ' ...
    'constant="false" complex="true">' ] );
writeParameter( fileID, 'name', model.objectModel.stat.name, 7, 'false', 'false', ...
    'string' );
writeParameter( fileID, 'beta', model.objectModel.stat.beta, 7, 'false', 'false', ...
    'double' );
fprintf( fileID, '\t\t\t\t\t\t\t%s\n', [ '<parameter name="transfrom" ' ...
    'constant="false" complex="true">' ] );
writeParameter( fileID, 'funname', model.objectModel.stat.transform.funname, ...
    8, 'false', 'false', 'string', 'function' );
fprintf( fileID, '\t\t\t\t\t\t\t%s\n', '</parameter>' );
fprintf( fileID, '\t\t\t\t\t\t%s\n', '</parameter>' );

fprintf( fileID, '\t\t\t\t\t\t%s\n', [ '<parameter name="relation" ' ...
    'constant="false" complex="true">' ] );
writeParameter( fileID, 'funname', model.objectModel.relation.funname, 7, ...
    'false', 'false', 'string', 'function' );
writeParameter( fileID, 'nvar', model.objectModel.relation.nvar, ...
    7, 'true', 'false', 'string' );
fprintf( fileID, '\t\t\t\t\t\t\t%s\n', [ '<parameter name="param" '...
    'constant="false" complex="false" type="none" structure="cell">' ]);
fprintf( fileID, '\t\t\t\t\t\t\t\t%s\n', [ '<parameter name="" ' ...
    'constant="false" complex="true" type="none" structure="struct">' ]);
writeParameter( fileID, 'funname', model.objectModel.relation.param{1}.funname, ...
    9, 'false', 'false', 'string', 'function' );
writeParameter( fileID, 'nvar', model.objectModel.relation.param{1}.nvar, ...
    9, 'true', 'false', 'double' );
fprintf( fileID, '\t\t\t\t\t\t\t\t%s\n', '</parameter>' );
fprintf( fileID, '\t\t\t\t\t\t\t\t%s\n', [ '<parameter name="" ' ...
    'constant="false" complex="true" type="none" structure="struct">' ]);
writeParameter( fileID, 'funname', model.objectModel.relation.param{2}.funname, ...
    9, 'false', 'false', 'double' );
fprintf( fileID, '\t\t\t\t\t\t\t\t%s\n', '</parameter>' );
fprintf( fileID, '\t\t\t\t\t\t\t%s\n', '</parameter>' );
fprintf( fileID, '\t\t\t\t\t\t%s\n', '</parameter>' );

fprintf( fileID, '\t\t\t\t\t\t%s\n', [ '<parameter name="intensStatModel" ' ...
    'constant="false" complex="true">' ] );
writeParameter( fileID, 'name', model.objectModel.intensStatModel.name, 7, ...
    'false', 'false', 'string' );
writeParameter( fileID, 'mu', model.objectModel.intensStatModel.mu, 7, ...
    'true', 'false', 'double' );
writeParameter( fileID, 'sigma', model.objectModel.intensStatModel.sigma, 7, ...
    'true', 'false', 'double' );
fprintf( fileID, '\t\t\t\t\t\t%s\n', '</parameter>' );

fprintf( fileID, '\t\t\t\t\t\t%s\n', [ '<parameter name="numStatModel" ' ...
    'constant="false" complex="true">' ] );
writeParameter( fileID, 'name', model.objectModel.numStatModel.name, 7, ...
    'false', 'false', 'string' );
writeParameter( fileID, 'alpha', model.objectModel.numStatModel.alpha, 7, ...
    'true', 'false', 'double' );
writeParameter( fileID, 'beta', model.objectModel.numStatModel.beta, 7, ...
    'true', 'false', 'double' );
fprintf( fileID, '\t\t\t\t\t\t%s\n', '</parameter>' );
fprintf( fileID, '\t\t\t\t\t%s\n', '</object>' );

%position model
fprintf( fileID, '\t\t\t\t\t%s\n', '<position>' );
writeParameter( fileID, 'name', model.positionModel.name, 6, 'false', ...
    'false', 'string' );
writeParameter( fileID, 'beta', model.positionModel.beta, 6, 'true', ...
    'false', 'double', 'array' );
fprintf( fileID, '\t\t\t\t\t\t%s\n', [ '<parameter name="transform" ' ...
    'constant="false" complex="false">' ] );
writeParameter( fileID, 'funname', model.positionModel.transform.funname, 7, ...
    'false', 'false', 'string', 'function' );
fprintf( fileID, '\t\t\t\t\t\t\t%s\n', [ '<parameter name="param" ' ...
    'constant="false" complex="false">'] );
writeParameter( fileID, 'order', model.positionModel.transform.param.order, ...
    8, 'true', 'false', 'double' );
writeParameter( fileID, 'scale', model.positionModel.transform.param.scale, ...
    8, 'true', 'false', 'double' );
fprintf( fileID, '\t\t\t\t\t\t\t%s\n', '</parameter>' );
fprintf( fileID, '\t\t\t\t\t\t%s\n', '</parameter>' );
fprintf( fileID, '\t\t\t\t\t%s\n', '</position>' );
fprintf( fileID, '\t\t\t\t%s\n', '</compartment>' );
end%proteinModel2slml

function proteinClass = getClass( compartmentName )
%GETCLASS Helper function used to determine the protein class as an input
%from the user at instance creation
switch lower( compartmentName )
    case{ '', 'unknown' }
        proteinClass = '0';
    case 'nucleus'
        proteinClass = '1';
    case{ 'endosome', 'endo' }
        proteinClass = '2';
    case{ 'golgi', 'golgi apparatus' }
        proteinClass = '3';
    case{ 'lysosome', 'lyso' }
        proteinClass = '4';
    case{ 'mitochondria', 'mitochondrion', 'mito' }
        proteinClass = '5';
    case{ 'nucleolus', 'nuc' }
        proteinClass = '6';
    otherwise
        proteinClass = '0';
end
end%getClass