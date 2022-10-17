function model = slml2slml( varargin )
% SLML2SLML Combines multiple generative model files into a single model file.
%
% List Of Input Arguments     Descriptions
% -----------------------     ------------
% files                       list of paths of models need be combined
% options                     Options structure
%
% The input argument options holds the valid parameters for these components.
% The shape of options is described below
%
% List Of Parameters        Descriptions
% ------------------        ------------
% output_filename           (optional) the file name of output model.
%                           Default is "model.mat".
% selection                 (mandatory) a matrix used to specify what submodels
%                           should be used from each file.
%
% Documentation (optional)
% ------------------------
% This is an optional structure with multiple elements that holds documentation about this model.
% If the documentation is not input, function will inherit documentation from first
% model in list if model is present.

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
%
% Copyright (C) 2018-2019 Murphy Lab
% Computational Biology Department
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

model = struct([]);
if isdeployed
    
    filename = is_deployed(varargin);
    load(filename);
else
    % @icaoberg method must accept exactly two parameters
    if nargin ~= 2
        warning('Wrong number of input arguments. Exiting method');
        return
    end

    files = varargin{1};
    options = varargin{2};

    if isempty( files )
        warning('Input argument files cannot be empty.')
        return
    end
end

if length( files ) > 3
    warning(['The currently implementation only allows the ' ...
        'concatenation of 3 generative model files. Ignoring models on index > 3']);
    files = files(1:3);
end

if isempty( options )
    warning('Options argument cannot be empty.');
    return
end

if ~isfield( options, 'output_filename' )
    options.output_filename = [ pwd filesep 'model.mat' ];
end

if length( files ) == 1
    load(files{1})
    save( options.output_filename, 'model' );
    return
end

if ~isfield( options, 'selection' )
    warning( 'Input argument selection cannot be empty. Assuming default values.' )
    if length(files) == 2
        options.selection = [1,1,0;0,0,1];
    end

    if length(files) == 3
        options.selection = eye(3);
    end

    if length(files) ~= 2 && length(files) ~= 3
        warning(['Length of input arguments is ' num2str(length(file)) ...
            '. Unable to assume default values.']);
    end
end

if ~isfield( options, 'documentation' )
    disp(['Documentation is not present in options structure.' ...
        ' Method will inherit documentation from first model in list if model is present.']);
end

for i=1:1:length(files)
    if ~exist( files{i} )
        warning(['File ' files{i} ' does not exist. Exiting method']);
        return
    end
end

models = {};
for i=1:1:length(files)
    models{i} = load( files{i} );
end

dimensionality = models{1}.model.dimensionality;
for i=2:1:length(models)
    if ~strcmpi( dimensionality, models{i}.model.dimensionality )
        warning('All models should have the same dimensionality. Exiting method.')
        return
    end
end

try
    temporary_model.nuclearShapeModel = ...
        models{find(options.selection(:,1)==1)}.model.nuclearShapeModel;
catch err
    getReport(err)
    warning('Unable to extract nuclear shape model. Exiting method.');
    return
end

try
    temporary_model.cellShapeModel = ...
        models{find(options.selection(:,2)==1)}.model.cellShapeModel;
catch err
    getReport(err)
    warning('Unable to extract cell shape model. Exiting method.');
    return
end

try
    temporary_model.proteinModel = ...
        models{find(options.selection(:,3)==1)}.model.proteinModel;
catch err
    getReport(err)
    warning('Unable to extract protein model. Exiting method.');
    return
end

if isfield( options, 'documentation' )
    temporary_model.documentation = options.documentation;
end

if ~isfield( options, 'documentation' ) && ...
        isfield( models{1}.model, 'documentation' )
    temporary_model.documentation = models{1}.model.documentation;
end

if isfield( options, 'name' )
    temporary_model.name = options.name;
end

temporary_model.dimensionality = models{1}.model.dimensionality;
temporary_model.id = uuidgen();
temporary_model.nuclearShapeModel.id = uuidgen();
temporary_model.cellShapeModel.id = uuidgen();
temporary_model.proteinModel.id = uuidgen();
temporary_model.documentation.date = date;
temporary_model.documentation.original_files = files;
temporary_model.documentation.selection = options.selection;

model = temporary_model; clear temporary_model;
save( options.output_filename, 'model' );
end%slml2slml
