function [ options ] = nameModels( CSGdata, meshData, models, imgs, options )
%NAMEMODELS Name models for use in convertCSGToMesh and readNetworkIntoGeometry.
%
% Inputs
% ------
% CSGdata  = 
% meshData = 
%
% Outputs
% -------
% options = options with field options.output.model_names


% Authors: Taraz Buck
%
% Copyright (C) 2019-2020 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
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


warning_id = ['CellOrganizer:', mfilename];
warning('off', warning_id);

if ~isfield(options.output, 'use_individual_meshes')
    options.output.use_individual_meshes = false;
    % options.output.use_individual_meshes = true;
end
use_individual_meshes = options.output.use_individual_meshes;


% Names of loaded models other than cell and nucleus

model_names = cell(length(models), 1);
model_cytonuclearflags = cell(length(models), 1);
for i = 1:length(models)
    model = models{i};
    if isfield(model, 'filename')
        % Assume no duplicated filenames, make insensitive to input order
        name = strrep(model.filename, '.', '_');
    elseif isfield(model.documentation, 'original_files')
        name = [];
        for j = 1:length(model.documentation.original_files)
            [original_file_path, original_file_name, original_file_ext] = fileparts(model.documentation.original_files{j});
            if j > 1
                name = [name, '_'];
            end
            name = [name, original_file_name, '_', strrep(original_file_ext, '.', '')];
        end
    else
        name = ['model', num2str(i)];
    end
    cytonuclearflag = 'all';
    if isfield(model, 'proteinModel') && isfield(model.proteinModel, 'cytonuclearflag')
        cytonuclearflag = model.proteinModel.cytonuclearflag;
    end
    model_names{i} = name;
    model_cytonuclearflags{i} = cytonuclearflag;
end
options.output.model_names = model_names;
options.output.model_cytonuclearflags = model_cytonuclearflags;


end
