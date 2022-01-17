function [ network_with_geometry_info ] = readNetworkIntoGeometry( CSGdata, meshData, models, imgs, options )
%READNETWORKINTOGEOMETRY Combines CellOrganizer geometry with compartmental BioNetGen-generated NET file reaction network.
%
% Inputs
% ------
% models   = cell array of models produced in slml2img.m
% imgs     = cell array containing image arrays for all objects to be saved
% options  = options struct produced in slml2img.m
%
% Outputs
% -------
% result = boolean flag indicating success
%
% Notes
% -----
% * options.output.NET should be the network file generated from a cBNG (compartmental BioNetGen) file by BNG2.pl.
% * "Bridged-surface connected transport" does not appear to be supported in VCell (see https://github.com/virtualcell/vcell/blob/088cd0a3164addefde730b046f269952056699d8/vcell-core/src/main/java/cbit/vcell/mapping/AbstractMathMapping.java#L1363).
% * Certain reactions between components in different compartments do not appear to be supported in VCell.
% * Tested only with a NET file with reactions in one compartment.


% Authors: Taraz Buck, Rohan Arepally, and Devin Sullivan
%
% Copyright (C) 2012-2020 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
% Copyright (C) 2018-2019 Taraz Buck
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

% 2019-06-11 Taraz Buck - Copied `instance2VCML.m` to `readNetworkIntoGeometry.m`.


warning('off', 'CellOrganizer:readNetworkIntoGeometry');

if nargin < 3
    warning('CellOrganizer:readNetworkIntoGeometry', 'First three arguments are required');
    return;
end

imgs_valid = true(length(imgs), 1);
for i = 1:length(imgs)
    imgs_valid(i) = imgs_valid(i) && numel(imgs{i}) > 0;
    imgs_valid(i) = imgs_valid(i) && ndims(imgs{i}) == 3;
end
if ~all(imgs_valid)
    warning('CellOrganizer:readNetworkIntoGeometry', 'Argument imgs must contain nonempty 3D arrays');
    return;
end


network_info = readNetwork(options.output.NET.filename);


network_info_parameters_keys = [];
parameters_names_to_expressions_map = [];
parameters_names_topological_order = [];

% Convert network_info.parameters to a map from name to DimensionedExpression
network_info_parameters_keys = network_info.parameters.keys();
for i = 1:length(network_info_parameters_keys)
    parameter_name = network_info_parameters_keys{i};
    parameter = network_info.parameters(parameter_name);
    parameter_value_expression = parameter.value_expression;
    parameter = DimensionedExpression(parameter_value_expression);
    network_info.parameters(parameter_name) = parameter;
end

parameters_should_write_to_parameters = containers.Map('KeyType', 'char', 'ValueType', 'uint8');
function set_should_write_to_parameters(given_parameter_name, given_parameter_should_write_to_parameters, given_parameters_should_write_to_parameters)
    if nargin < 3
        given_parameters_should_write_to_parameters = parameters_should_write_to_parameters;
    end
    if ~given_parameter_should_write_to_parameters
        given_parameters_should_write_to_parameters(given_parameter_name) = false;
    else
        if given_parameters_should_write_to_parameters.isKey(given_parameter_name)
            given_parameters_should_write_to_parameters.remove(given_parameter_name);
        end
    end
end
function [given_parameter_should_write_to_parameters] = get_should_write_to_parameters(given_parameter_name, given_parameters_should_write_to_parameters)
    if nargin < 2
        given_parameters_should_write_to_parameters = parameters_should_write_to_parameters;
    end
    if given_parameters_should_write_to_parameters.isKey(given_parameter_name) && given_parameters_should_write_to_parameters(given_parameter_name) == false
        given_parameter_should_write_to_parameters = false;
    else
        given_parameter_should_write_to_parameters = true;
    end
end

function updateParametersMap(given_parameter_name, given_parameter, given_parameter_should_write_to_parameters)
    if nargin < 3
        given_parameter_should_write_to_parameters = true;
    end
    if nargin >= 2
        if ~network_info.parameters.isKey(given_parameter_name)
            network_info_parameters_keys{end+1} = network_info_parameters_keys;
            % parameters_names_topological_order{end+1} = given_parameter_name;
            parameters_names_topological_order = [];
        end
        network_info.parameters(given_parameter_name) = given_parameter;
        set_should_write_to_parameters(given_parameter_name, given_parameter_should_write_to_parameters);
    %{
    else
        network_info_parameters_keys = network_info.parameters.keys();
        % parameters_names_to_expressions_map = containers.Map(network_info.parameters.keys(), cellfun(@(x)x.value, network_info.parameters.values(), 'UniformOutput', false));
        parameters_names_to_expressions_map = containers.Map(network_info.parameters.keys(), network_info.parameters.values());
    %}
    else
        network_info_parameters_keys = network_info.parameters.keys();
        parameters_names_topological_order = topologicalSort(network_info.parameters, @(x)network_info.parameters(x).getVariables());
    end
    parameters_names_to_expressions_map = network_info.parameters;
end

function [varargout] = inferUnits(given_expression)
    % Determine units to the extent possible
    if nargin >= 1
        % if ~ischar(given_expression)
        if isa(given_expression, 'DimensionedExpression')
            given_expression = given_expression.value;
        end
        if ~ischar(given_expression)
            error('given_expression must be string or DimensionedExpression');
        end
        given_expression = char(given_expression);
        given_expression_evaluated = mathEval(given_expression, parameters_names_to_expressions_map, struct('check_variables', false, 'return_type', 'DimensionedExpression'));
        result_expression = DimensionedExpression(given_expression, given_expression_evaluated.units_manager);
        varargout = {result_expression};
    else
        updateParametersMap();
        for i = 1:length(network_info_parameters_keys)
            % parameter_name = network_info_parameters_keys{i};
            parameter_name = parameters_names_topological_order{i};
            parameter = network_info.parameters(parameter_name);
            % if parameter.isDimensionless()
            if parameter.hasVariables()
                parameter = inferUnits(parameter);
                updateParametersMap(parameter_name, parameter);
            end
        end
        varargout = {};
    end
end

function result = parameterReferenceBase(given_parameter_name, given_parameters)
    result = DimensionedExpression(given_parameter_name, given_parameters(given_parameter_name).units_manager);
end

function result = parameterReference(given_parameter_name)
    result = parameterReferenceBase(given_parameter_name, network_info.parameters);
end

function result = evaluateExpression(given_expression, given_parameters_names_to_expressions_map)
    if nargin < 2
        given_parameters_names_to_expressions_map = parameters_names_to_expressions_map;
    end
    given_expression = DimensionedExpression(given_expression);
    result = mathEval(given_expression, given_parameters_names_to_expressions_map, struct('check_variables', true, 'return_type', 'DimensionedExpression'));
    result = DimensionedExpression(result.value, given_expression.units_manager);
end

function printParameterValues()
    fprintf('\n')
    fprintf('*** Parameter values (%d parameters) ***\n', length(network_info_parameters_keys))
    for i = 1:length(network_info_parameters_keys)
        parameter_name = network_info_parameters_keys{i};
        parameter = network_info.parameters(parameter_name);
        parameter_evaluated = evaluateExpression(parameter);
        fprintf('*** %s = %s = %s\n', parameter_name, char(parameter), char(parameter_evaluated))
    end
    fprintf('\n')
end

% updateParametersMap();
inferUnits();



warning('CellOrganizer:readNetworkIntoGeometry', 'Change options naming for these variables!');
translations = options.output.NET.translations;
model_names = options.output.model_names;
model_cytonuclearflags = options.output.model_cytonuclearflags;
downsampling = options.output.NET.downsampling;
if length(downsampling) == 1
    downsampling = downsampling .* [1, 1, 1];
end
%{
if downsampling(1) ~= downsampling(2) || downsampling(1) ~= downsampling(3)
    error('downsampling is not uniform! Implement and use DimensionedArray if required.');
end
downsampling_single = downsampling(1);
%}
add_translocation_intermediates = options.output.NET.add_translocation_intermediates;

default_diffusion_coefficient = options.output.NET.default_diffusion_coefficient;

net_length_unit = UnitsManager(options.output.NET.units.length);
net_area_unit = net_length_unit^2;
net_volume_unit = net_length_unit^3;
net_time_unit = UnitsManager(options.output.NET.units.time);
net_count_unit = UnitsManager('molecule');
net_concentration_unit = UnitsManager(options.output.NET.units.concentration);
net_diffusion_coefficient_unit = net_area_unit / net_time_unit;
warning('CellOrganizer:readNetworkIntoGeometry', 'Assuming NET file concentrations are extensive');

avogadro_constant = 6.022140857e23;


output_NET_defaults.output_length_unit = 'um';
output_NET_defaults.output_concentration_unit = 'uM';
output_NET_defaults.output_membrane_substance_unit = 'molecule';
output_NET_defaults.output_time_unit = 's';
output_NET_defaults.use_image_adjacency = false;
% output_NET_defaults.use_image_adjacency = true;
if ~isfield( options.output, 'NET' )
    options.output.NET = output_NET_defaults;
else
    options.output.NET = ml_initparam( options.output.NET, output_NET_defaults );
end

% use_image_adjacency = options.output.NET.use_image_adjacency;

output_length_unit = UnitsManager(options.output.NET.output_length_unit);
output_concentration_unit = UnitsManager(options.output.NET.output_concentration_unit);
output_membrane_substance_unit = UnitsManager(options.output.NET.output_membrane_substance_unit);
output_time_unit = UnitsManager(options.output.NET.output_time_unit);

output_area_unit = output_length_unit^2;
output_volume_unit = output_length_unit^3;
% output_volume_substance_unit = output_concentration_unit * output_volume_unit;
output_volume_reaction_rate_unit = output_concentration_unit / output_time_unit;
% output_membrane_reaction_rate_unit = output_membrane_substance_unit / output_time_unit^2 / output_time_unit;
output_diffusion_coefficient_unit = output_area_unit / output_time_unit;

warning('CellOrganizer:readNetworkIntoGeometry', 'Assuming options.resolution is in um');
if options.resolution.cubic(1) ~= options.resolution.cubic(2) || options.resolution.cubic(1) ~= options.resolution.cubic(3)
    error('options.resolution.cubic is not uniform! Implement and use DimensionedArray if required.');
end
resolution = double(UnitsManager('um') / output_length_unit) * options.resolution.cubic;
% resolution_single = unit_convert('um', output_length_unit, resolution(1));
% resolution_single = double(output_length_unit / UnitsManager('um')) * resolution(1);
resolution_single = resolution(1);
voxel_volume = resolution_single^3;
voxel_face_area = resolution_single^2;

SI_effective_width = DimensionedExpression(options.output.NET.effective_width, 'm');
output_effective_width = SI_effective_width;

% use_image_adjacency = options.output.NET.useImageAdjacency;
use_image_adjacency = false;
% use_image_adjacency = true;

is_named_imgs_inclusive = true;

enforce_parent_buffer = true;
enforce_sibling_buffer = true;



% Constants

dim_chars = 'XYZ';
sign_chars = 'mp';
connectivity = 6;
if connectivity == 6
    connectivity_se = cat(3, [0, 0, 0; 0, 1, 0; 0, 0, 0], [0, 1, 0; 1, 1, 1; 0, 1, 0], [0, 0, 0; 0, 1, 0; 0, 0, 0]);
elseif connectivity == 18
    connectivity_se = cat(3, [0, 1, 0; 1, 1, 1; 0, 1, 0], [1, 1, 1; 1, 1, 1; 1, 1, 1], [0, 1, 0; 1, 1, 1; 0, 1, 0]);
elseif connectivity == 26
    connectivity_se = cat(3, [1, 1, 1; 1, 1, 1; 1, 1, 1], [1, 1, 1; 1, 1, 1; 1, 1, 1], [1, 1, 1; 1, 1, 1; 1, 1, 1]);
end

compartment_volume_prefix = 'vol_';
compartment_area_prefix = 'sa_';
compartment_diffusion_coefficient_prefix = 'dc_';
backup_prefix = 'backup_';

% kinetics_type = 'GeneralKinetics';
kinetics_type = 'MassAction';

avogadro_constant_value_expression = '6.022140857e23';
avogadro_constant_units = 'molecule.mol-1';

% output_diffusion_coefficient = unit_convert('m2.s-1', output_diffusion_coefficient_unit, default_diffusion_coefficient);
% output_diffusion_coefficient = DimensionedExpression(default_diffusion_coefficient, 'm2.s-1');
output_diffusion_coefficient = DimensionedExpression(default_diffusion_coefficient, net_diffusion_coefficient_unit);




% Names of loaded models other than cell and nucleus

function x = isStringNumeric(x)
    x = ~isnan(str2double(x));
end


function [given_pair_string_translated, given_pair_string_cell, given_pair_string] = get_compartment_pair_string(given_name1, given_name2, given_name_map)
    given_pair_string_cell = {translateWithDefaultIdentity(given_name_map, given_name1), translateWithDefaultIdentity(given_name_map, given_name2)};
    given_pair_string_cell = sort(given_pair_string_cell);
    given_pair_string = [given_pair_string_cell{1}, '_', given_pair_string_cell{2}];
    given_pair_string_translated = translateWithDefaultIdentity(given_name_map, given_pair_string);
end


function [given_named_imgs, given_named_imgs_exclusive_volumes] = get_exclusive_images(given_all_compartments_image, given_voxel_volume, given_names_sorted)
    % Separate indexed regions
    given_named_imgs = containers.Map('KeyType', 'char', 'ValueType', 'any');
    given_named_imgs_exclusive_volumes = containers.Map('KeyType', 'char', 'ValueType', 'double');
    for gi = 1:length(given_names_sorted)
        given_name = given_names_sorted{gi};
        given_img = given_all_compartments_image == gi;
        given_named_imgs(given_name) = given_img;
        given_named_imgs_exclusive_volumes(given_name) = sum(given_img(:)) * given_voxel_volume;
    end
end


function [given_named_imgs, given_named_imgs_exclusive_volumes, given_names_sorted] = separate_compartment_images(given_all_compartments_image, given_voxel_volume, given_names_sorted, given_network_info_compartments)
    % Recreate named_imgs by separating indexed regions in given_all_compartments_image and filling blank portions of parents with child regions
    
    [given_named_imgs, given_named_imgs_exclusive_volumes] = get_exclusive_images(given_all_compartments_image, given_voxel_volume, given_names_sorted);
    
    % Fill holes in parents using children (named_imgs are inclusive)
    for gi = length(given_names_sorted):-1:1
        given_name = given_names_sorted{gi};
        given_img = given_named_imgs(given_name);
        given_outside_name = given_network_info_compartments(given_name).outside;
        if ~isempty(given_outside_name)
            given_outside_name = given_network_info_compartments(given_outside_name).outside;
            given_named_imgs(given_outside_name) = given_named_imgs(given_outside_name) | given_img;
        end
    end
end


function [given_all_compartments_image, given_named_imgs, given_names_sorted] = combine_compartment_images(given_named_imgs, given_names_sorted, given_voxel_volume, given_network_info_compartments)
    % Combine compartments of named_imgs into an indexed image
    
    % Combine compartments with optional buffer from parent compartment boundary
    given_all_compartments_image = [];
    gk = 1;
    given_names_sorted2 = {};
    for gi = 1:length(given_named_imgs)
        given_name = given_names_sorted{gi};
        given_img = given_named_imgs(given_name);
        
        % Produce compartment mask
        given_mask = given_img;
        
        % Remove voxels adjacent to outside compartments' boundaries
        given_compartment = given_network_info_compartments(given_name);
        if ~isempty(given_compartment.outside)
            given_compartment_outside_name = given_network_info_compartments(given_compartment.outside).outside;
            
            given_compartment_outside_mask = given_named_imgs(given_compartment_outside_name);
            if enforce_parent_buffer
                given_compartment_outside_mask = imerode(given_compartment_outside_mask, connectivity_se);
            end
            given_mask(~given_compartment_outside_mask) = false;
        end
        
        if isempty(given_all_compartments_image)
            given_all_compartments_image = zeros(size(given_mask), 'uint8');
        end
        
        % Overwrite lower priority compartments like EC and framework
        % EC will be one, cell with be 2, etc.
        given_mask_mean = mean(given_mask(:));
        if given_mask_mean > 0
            given_all_compartments_image(given_mask) = gk;
            gk = gk + 1;
            given_names_sorted2{end+1} = given_name;
        end
        
        if given_mask_mean == 0
            % error('given_mask_mean == 0');
            warning('CellOrganizer:instance2VCML', 'Compartment ''%s'' given_mask_mean == 0', given_name);
        end
    end
    
    given_names_sorted = given_names_sorted2;
    
    [given_named_imgs, given_named_imgs_exclusive_volumes, given_names_sorted] = separate_compartment_images(given_all_compartments_image, given_voxel_volume, given_names_sorted, given_network_info_compartments);
    
    
    if enforce_sibling_buffer
        % Optional buffer from sibling compartment boundaries
        given_all_compartments_image2 = zeros(size(given_mask), 'uint8');
        gk = 1;
        given_names_sorted2 = {};
        for gi = 1:length(given_named_imgs)
            given_name = given_names_sorted{gi};
            given_img = given_named_imgs(given_name);
            
            % Produce compartment mask
            given_mask = given_img;
            
            % Remove voxels adjacent to outside compartments' boundaries
            given_compartment = given_network_info_compartments(given_name);
            if ~isempty(given_compartment.outside)
                % Get volume outside surrounding membrane
                given_compartment_outside_name = given_network_info_compartments(given_compartment.outside).outside;
                % Get membranes within outside volume
                given_compartment_sibling_names = given_network_info_compartments(given_compartment_outside_name).insides;
                % Get volumes within membranes
                for gj = 1:length(given_compartment_sibling_names)
                    given_compartment_sibling_names{gj} = given_network_info_compartments(given_compartment_sibling_names{gj}).insides{1};
                end
                % Remove self
                given_compartment_sibling_names = given_compartment_sibling_names(~strcmp(given_compartment_sibling_names, given_name));
                
                if ~isempty(given_compartment_sibling_names)
                    given_compartment_siblings_mask = false(size(given_mask));
                    for gj = 1:length(given_compartment_sibling_names)
                        given_compartment_sibling_mask = given_named_imgs(given_compartment_sibling_names{gj});
                        given_compartment_siblings_mask = given_compartment_siblings_mask | given_compartment_sibling_mask;
                    end
                    given_compartment_siblings_mask = imdilate(given_compartment_siblings_mask, connectivity_se);
                    given_mask(given_compartment_siblings_mask) = false;
                end
            end
            
            if isempty(given_all_compartments_image)
                given_all_compartments_image = zeros(size(given_mask), 'uint8');
            end
            
            % Overwrite lower priority compartments like EC and framework
            % EC will be one, cell with be 2, etc.
            given_mask_mean = mean(given_mask(:));
            if given_mask_mean > 0
                given_all_compartments_image2(given_mask) = gk;
                gk = gk + 1;
                given_names_sorted2{end+1} = given_name;
            end
            
            if given_mask_mean == 0
                % error('given_mask_mean == 0');
                warning('CellOrganizer:instance2VCML', 'Compartment ''%s'' given_mask_mean == 0', given_name);
            end
        end
        
        given_names_sorted = given_names_sorted2;
        given_all_compartments_image = given_all_compartments_image2;
        
        [given_named_imgs, given_named_imgs_exclusive_volumes, given_names_sorted] = separate_compartment_images(given_all_compartments_image, given_voxel_volume, given_names_sorted, given_network_info_compartments);
    end
    
    
    given_all_compartments_image_zero_mean = mean(given_all_compartments_image(:) == 0);
end


function [given_image_boundary_value] = get_compartment_image_boundary_value(given_image)
    given_image_boundary_values = {};
    given_image_boundary_values{end+1} = given_image(:, :, 1);
    given_image_boundary_values{end+1} = given_image(:, :, end);
    given_image_boundary_values{end+1} = given_image(:, 1, 2:end-1);
    given_image_boundary_values{end+1} = given_image(:, end, 2:end-1);
    given_image_boundary_values{end+1} = given_image(1, 2:end-1, 2:end-1);
    given_image_boundary_values{end+1} = given_image(end, 2:end-1, 2:end-1);
    given_image_boundary_values = cellfun(@(x)reshape(x, 1, []), given_image_boundary_values, 'UniformOutput', false);
    given_image_boundary_values = cell2mat(given_image_boundary_values);
    given_image_boundary_values_unique = unique(given_image_boundary_values);
    if length(given_image_boundary_values_unique) > 1
        error('Boundary is not uniform');
    end
    given_image_boundary_value = given_image_boundary_values_unique(1);
end


function [given_image] = crop_compartment_image(given_image)
    given_image_boundary_value = get_compartment_image_boundary_value(given_image);
    given_image_not_boundary = given_image ~= given_image_boundary_value;
    [~, given_image_crop_bounds] = cropImg(given_image_not_boundary, 1);
    given_image = given_image(given_image_crop_bounds(1):given_image_crop_bounds(2), given_image_crop_bounds(3):given_image_crop_bounds(4), given_image_crop_bounds(5):given_image_crop_bounds(6));
end


function result = combine_resize_compartment_images(given_named_imgs, given_voxel_volume, given_names_sorted, given_network_info_compartments, given_scale)
    % Combine named_imgs into an indexed image, resample, and reproduce named_imgs
    
    [given_all_compartments_image, given_named_imgs, given_names_sorted] = combine_compartment_images(given_named_imgs, given_names_sorted, given_voxel_volume, given_network_info_compartments);
    
    given_scale_pad_size = ceil(1 ./ given_scale);
    
    given_all_compartments_image = padarray(given_all_compartments_image, given_scale_pad_size, 'replicate');
    given_all_compartments_image_before_downsampling_size = size(given_all_compartments_image);
    given_all_compartments_image = image_resize_nd(given_all_compartments_image, given_scale, 'nearest');
    given_all_compartments_image_size = size(given_all_compartments_image);
    
    given_actual_scale = given_all_compartments_image_size ./ given_all_compartments_image_before_downsampling_size;
    given_voxel_volume = given_voxel_volume ./ prod(given_actual_scale);
    
    % Assumes boundary should be one value and is of arbitrary size
    given_all_compartments_image = crop_compartment_image(given_all_compartments_image);
    
    [given_named_imgs, given_named_imgs_exclusive_volumes, given_names_sorted] = separate_compartment_images(given_all_compartments_image, given_voxel_volume, given_names_sorted, given_network_info_compartments);
    [given_all_compartments_image, given_named_imgs, given_names_sorted] = combine_compartment_images(given_named_imgs, given_names_sorted, given_voxel_volume, given_network_info_compartments);
    [given_named_imgs, given_named_imgs_exclusive_volumes, given_names_sorted] = separate_compartment_images(given_all_compartments_image, given_voxel_volume, given_names_sorted, given_network_info_compartments);
    
    result = struct();
    result.all_compartments_image = given_all_compartments_image;
    result.named_imgs = given_named_imgs;
    result.actual_scale = given_actual_scale;
    result.named_imgs_exclusive_volumes = given_named_imgs_exclusive_volumes;
    result.names_sorted = given_names_sorted;
end


% Translate names

% Create a map for translating compartment names
name_map = containers.Map('KeyType', 'char', 'ValueType', 'char');
name_map('EC') = 'EC';
name_map('cell') = 'cell';
name_map('nucleus') = 'nucleus';
name_map('nuc') = 'nucleus';
% Defaults
for i = 1:length(models)
    name_map(model_names{i}) = model_names{i};
end
% Translations in options
for i = 1:size(translations, 1)
    name_map(translations{i, 1}) = translations{i, 2};
end

EC_translated = translateWithDefaultIdentity(name_map, 'EC');
cell_translated = translateWithDefaultIdentity(name_map, 'cell');
nucleus_translated = translateWithDefaultIdentity(name_map, 'nucleus');
EC_cell_translated = get_compartment_pair_string('EC', 'cell', name_map);
cell_nucleus_translated = get_compartment_pair_string('cell', 'nucleus', name_map);
framework_compartment_names = {EC_translated, cell_translated, nucleus_translated};
is_framework_compartment_name_function = @(x)any(strcmp(x, framework_compartment_names));
% model_names = cellfun(@(x)translateWithDefaultIdentity(name_map, x), model_names, 'UniformOutput', false);

model_names_translated = cell(length(models), 1);
for i = 1:length(model_names)
    model_names_translated{i} = translateWithDefaultIdentity(name_map, model_names{i});
end

% Rename compartments and their fields
network_info_compartments_keys = network_info.compartments.keys;
for i = 1:length(network_info_compartments_keys)
    compartment_name = network_info_compartments_keys{i};
    compartment = network_info.compartments(compartment_name);
    compartment_name_translated = translateWithDefaultIdentity(name_map, compartment_name);
    compartment.name = compartment_name_translated;
    compartment.outsides = translateWithDefaultIdentity(name_map, compartment.outsides);
    for j = 1:length(compartment.insides)
        compartment.insides{j} = translateWithDefaultIdentity(name_map, compartment.insides{j});
    end
    network_info.compartments.remove(compartment_name);
    network_info.compartments(compartment_name_translated) = compartment;
end
network_info_compartments_keys = network_info.compartments.keys;

% Rename meshes
for i = 1:length(meshData)
    meshData(i).source = 'mesh';
    for j = 1:length(meshData(i).list)
        meshData_list_item_name = meshData(i).list(j).name;
        meshData_list_item_name = translateWithDefaultIdentity(name_map, meshData_list_item_name);
        meshData(i).list(j).name = meshData_list_item_name;
    end
end

% Rename compartments
for i = 1:length(network_info.compartments)
    compartment_name = network_info_compartments_keys{i};
    network_info_compartments_keys
    meshData(i).source = 'mesh';
    for j = 1:length(meshData(i).list)
        meshData_list_item_name = meshData(i).list(j).name;
        meshData_list_item_name = translateWithDefaultIdentity(name_map, meshData_list_item_name);
        meshData(i).list(j).name = meshData_list_item_name;
        network_info.compartments(meshData_list_item_name).outside = translateWithDefaultIdentity(name_map, network_info.compartments(meshData_list_item_name).outside);
        for k = 1:length(network_info.compartments(meshData_list_item_name).insides)
            network_info.compartments(meshData_list_item_name).insides{k} = translateWithDefaultIdentity(name_map, network_info.compartments(meshData_list_item_name).insides{k});
        end
    end
end

%{
% Rename meshes named for membranes after the volumes they contain
for i = 1:length(meshData)
    meshData(i).source = 'mesh';
    % meshData(i).list = setfield(meshData(i).list, '')
    meshData(i).list.surface_name = meshData(i).list.name;
    for j = 1:length(meshData(i).list)
        meshData_list_item_name = meshData(i).list(j).surface_name;
        if network_info.compartments.isKey(meshData_list_item_name) && network_info.compartments(meshData_list_item_name).spatial_dimensions == 2
            meshData_list_item_name = network_info.compartments(meshData_list_item_name).insides{1};
            meshData(i).list(j).surface_name = meshData_list_item_name;
        end
    end
end
%}



% Create a single image containing all compartments. Assumes no overlapping.

% Collect names and properties of compartments
warning('CellOrganizer:readNetworkIntoGeometry', 'Volumes might be handled differently than in BioNetGen and other software; compare https://github.com/RuleWorld/bionetgen/blob/master/bng2/Perl2/Compartment.pm');
named_imgs = containers.Map('KeyType', 'char', 'ValueType', 'any');
% Exclusive volumes (do not include volumes of compartments inside)
named_imgs_volumes = containers.Map('KeyType', 'char', 'ValueType', 'double');
named_imgs_cytonuclearflags = containers.Map('KeyType', 'char', 'ValueType', 'char');
% Topological sort, roughly in descending order by volume
names_sorted = {};
if any(strcmpi(options.synthesis, {'all', 'framework'}))
    img = imgs{2};
    if ~is_named_imgs_inclusive
        img = img & ~imgs{1};
    end
    named_imgs(cell_translated) = img;
elseif strcmpi(options.synthesis, 'cell')
    named_imgs(cell_translated) = imgs{1};
end
if any(strcmpi(options.synthesis, {'all', 'framework', 'cell'}))
    names_sorted{end+1} = cell_translated;
end
if any(strcmpi(options.synthesis, {'all', 'framework', 'nucleus'}))
    named_imgs(nucleus_translated) = imgs{1};
    names_sorted{end+1} = nucleus_translated;
end
if any(strcmpi(options.synthesis, {'all', 'cell'}))
    for i = 1:length(models)
        model = models{i};
        name = model_names_translated{i};
        cytonuclearflag = model_cytonuclearflags{i};
        named_imgs(name) = imgs{length(named_imgs)+1};
        named_imgs_cytonuclearflags(name) = cytonuclearflag;
        names_sorted{end+1} = name;
    end
end
names = named_imgs.keys;

% Make images logical
for i = 1:length(named_imgs)
    name = names{i};
    named_imgs(name) = named_imgs(name) > 0;
end

% Add padding where any shape touches the image boundary
padding_needed_before = zeros(1, 3);
padding_needed_after = zeros(1, 3);
for i = 1:length(named_imgs)
    name = names{i};
    img = named_imgs(name);
    padding_needed_before = max(padding_needed_before, [any(any(img(1, :, :))), any(any(img(:, 1, :))), any(any(img(:, :, 1)))]);
    padding_needed_after = max(padding_needed_before, [any(any(img(end, :, :))), any(any(img(:, end, :))), any(any(img(:, :, end)))]);
end
% Pad all images with empty space
padding_needed_before = padding_needed_before + 1;
padding_needed_after = padding_needed_after + 1;
for i = 1:length(named_imgs)
    name = names{i};
    img = named_imgs(name);
    img = padarray(img, padding_needed_before, 'pre');
    img = padarray(img, padding_needed_after, 'post');
    named_imgs(name) = img;
end

% Add extracellular space
ec_img = [];
for i = 1:length(named_imgs)
    name = names{i};
    named_img = named_imgs(name);
    if isempty(ec_img)
        ec_img = true(size(named_img));
    end
    if is_named_imgs_inclusive
        break;
    end
    ec_img(named_img) = false;
end
named_imgs(EC_translated) = ec_img;
names = named_imgs.keys;
names_sorted = [{EC_translated}, names_sorted];

% Collect names and properties of compartments

% Add missing framework compartments to network_info.compartments
if any(strcmpi(options.synthesis, {'all', 'framework', 'cell'}))
    if ~network_info.compartments.isKey(EC_translated)
        network_info.compartments(EC_translated) = struct('name', EC_translated, 'spatial_dimensions', 3, 'size_expression', nan, 'outside', '', 'insides', {{EC_cell_translated}}, 'double_outside', '', 'double_insides', {{cell_translated}});
    end
    if ~network_info.compartments.isKey(EC_cell_translated)
        network_info.compartments(EC_cell_translated) = struct('name', EC_cell_translated, 'spatial_dimensions', 2, 'size_expression', nan, 'outside', EC_translated, 'insides', {{cell_translated}}, 'double_outside', '', 'double_insides', {{}});
        if any(strcmpi(options.synthesis, {'all', 'framework'}))
            temp.double_insides = {cell_nucleus_translated};
        end
    end
    if ~network_info.compartments.isKey(cell_translated)
        temp = struct('name', cell_translated, 'spatial_dimensions', 3, 'size_expression', nan, 'outside', EC_cell_translated, 'insides', {{}}, 'double_outside', EC_translated, 'double_insides', {{}});
        if any(strcmpi(options.synthesis, {'all', 'framework'}))
            temp.insides = {cell_nucleus_translated};
            temp.double_insides = {nucleus_translated};
        end
        network_info.compartments(cell_translated) = temp;
    end
end
if ~network_info.compartments.isKey(nucleus_translated) && any(strcmpi(options.synthesis, {'all', 'framework', 'nucleus'}))
    temp = struct('name', nucleus_translated, 'spatial_dimensions', 3, 'size_expression', nan, 'outside', '', 'insides', {{}}, 'double_outside', '', 'double_insides', {{}});
    if any(strcmpi(options.synthesis, {'all', 'framework'}))
        temp.outside = cell_nucleus_translated;
        temp.double_outside = cell_translated;
    end
    network_info.compartments(nucleus_translated) = temp;
end
if ~network_info.compartments.isKey(cell_nucleus_translated) && any(strcmpi(options.synthesis, {'all', 'framework'}))
    network_info.compartments(cell_nucleus_translated) = struct('name', cell_nucleus_translated, 'spatial_dimensions', 2, 'size_expression', nan, 'outside', cell_translated, 'insides', {{nucleus_translated}}, 'double_outside', EC_cell_translated, 'double_insides', {{}});
end

% Add missing model compartments to network_info.compartments
if strcmpi(options.synthesis, 'all')
    for i = 1:length(models)
        model = models{i};
        name = model_names_translated{i};
        cytonuclearflag = model_cytonuclearflags{i};
        outside = '';
        outside_outside = '';
        cytonuclearflag = named_imgs_cytonuclearflags(name);
        switch cytonuclearflag
            case 'cyto'
                outside_outside = cell_translated;
                outside = get_compartment_pair_string(outside_outside, name, name_map);
            case 'nucleus'
                outside_outside = nucleus_translated;
                outside = get_compartment_pair_string(outside_outside, name, name_map);
            case 'all'
                error('Unable to assign outside to proteins/objects with cytonuclearflag ''all''');
        end
        
        if ~network_info.compartments.isKey(name)
            network_info.compartments(name) = struct('name', name, 'spatial_dimensions', 3, 'size_expression', nan, 'outside', outside, 'insides', {{}}, 'double_outside', '', 'double_insides', {{}});
        end
        
        if ~network_info.compartments.isKey(outside)
            network_info.compartments(outside) = struct('name', outside, 'spatial_dimensions', 2, 'size_expression', nan, 'outside', outside_outside, 'insides', {{name}}, 'double_outside', '', 'double_insides', {{}});
        end
        
        compartment_outside = network_info.compartments(outside);
        if ~any(strcmp(compartment_outside.insides, name))
            compartment_outside.insides{end+1} = name;
        end
        network_info.compartments(outside) = compartment_outside;
        
        compartment_outside_outside = network_info.compartments(outside_outside);
        if ~any(strcmp(compartment_outside_outside.double_insides, name))
            compartment_outside_outside.double_insides{end+1} = name;
        end
        network_info.compartments(outside_outside) = compartment_outside_outside;
    end
end


% Create combined image of all compartments. Each voxel will belong to only one compartment.
[all_compartments_image, named_imgs, names_sorted] = combine_compartment_images(named_imgs, names_sorted, voxel_volume, network_info.compartments);

for i = 1:length(names)
    name = names{i};
    img = named_imgs(name);
    % named_imgs_volumes(name) = sum(img(:)) * prod(resolution);
    named_imgs_volumes(name) = sum(img(:)) * voxel_volume;
end


n_compartments = length(names_sorted);
all_compartment_indices = containers.Map('KeyType', 'char', 'ValueType', 'double');
for i = 1:n_compartments
    all_compartment_indices(names_sorted{i}) = i;
end

% Add insides field to network_info.compartments values
n_network_info_compartments = length(network_info.compartments);
network_info_compartments_keys = network_info.compartments.keys;

for i = 1:n_network_info_compartments
    name = network_info_compartments_keys{i};
    compartment = network_info.compartments(name);
    if ~isfield(compartment, 'insides')
        compartment.insides = {};
    end
    network_info.compartments(name) = compartment;
end
for i = 1:n_network_info_compartments
    name = network_info_compartments_keys{i};
    outside_name = network_info.compartments(name).outside;
    if ~network_info.compartments.isKey(outside_name)
        continue;
    end
    outside_compartment = network_info.compartments(outside_name);
    if ~isempty(outside_name)
        outside_compartment.insides{end+1} = name;
    end
    network_info.compartments(outside_name) = outside_compartment;
end
for i = 1:n_network_info_compartments
    name = network_info_compartments_keys{i};
    compartment = network_info.compartments(name);
    if ~isfield(compartment, 'insides')
        continue;
    end
    compartment.insides = unique(compartment.insides);
    network_info.compartments(name) = compartment;
end


%{
% Create indexed image
all_compartments_image_before_downsampling = [];
all_compartments_image = [];
named_imgs_before_downsampling = containers.Map(named_imgs.keys(), named_imgs.values());
% resolution_before_downsampling = resolution;
resolution_single_before_downsampling = resolution_single;
voxel_volume_before_downsampling = voxel_volume;
% if any(downsampling ~= 1)
if downsampling_single ~= 1
    for i = 1:length(named_imgs)
        name = named_imgs_volumes_keys_sorted{i};
        img = named_imgs(name);
        % named_imgs_size_before = size(img)
        img = img > 0;
        % img = image_resize_nd(single(img), ones(1, 3) .* downsampling, 'bilinear') >= 0.5;
        img = image_resize_nd(single(img), ones(1, 3) * downsampling_single, 'bilinear') >= 0.5;
        img = padarray(img, ones(1, 3), false);
        named_imgs(name) = img;
        % named_imgs_size_after = size(img)
    end
    % resolution = resolution ./ downsampling;
    resolution_single = resolution_single / downsampling_single;
    voxel_volume = resolution_single / downsampling_single;
end
% options.output.resolution = resolution;
options.output.resolution = resolution_single .* [1, 1, 1];
for i = 1:length(named_imgs)
    name = named_imgs_volumes_keys_sorted{i};
    img_before_downsampling = named_imgs_before_downsampling(name);
    img = named_imgs(name);
    
    % Before downsampling
    mask = img_before_downsampling > 0;
    if isempty(all_compartments_image_before_downsampling)
        all_compartments_image_before_downsampling = zeros(size(mask), 'uint8');
    end
    % Don't overwrite objects with smaller total volumes
    mask = mask & (all_compartments_image_before_downsampling == 0);
    % EC will be zero, cell with be 1, etc.
    all_compartments_image_before_downsampling(mask) = i;
    mask_mean = mean(mask(:));
    
    % After downsampling
    mask = img > 0;
    if isempty(all_compartments_image)
        all_compartments_image = zeros(size(mask), 'uint8');
    end
    % Don't overwrite objects with smaller total volumes
    mask = mask & (all_compartments_image == 0);
    % EC will be zero, cell with be 1, etc.
    all_compartments_image(mask) = i;
    mask_mean = mean(mask(:));
    if mask_mean == 0
        % error('mask_mean == 0');
        warning('CellOrganizer:readNetworkIntoGeometry', 'Compartment ''%s'' mask_mean == 0', name);
    end
end


% Add extracellular space
n_compartments = length(names_sorted);
% Exclusive volumes (do not include volumes of compartments inside)
all_compartment_volumes = zeros(n_compartments, 1);
for i = 2:n_compartments
    all_compartment_volumes(i) = named_imgs_volumes(names_sorted{i});
end
named_imgs(EC_translated) = all_compartments_image == 0;
img = named_imgs(EC_translated);
% named_imgs_volumes(EC_translated) = sum(img(:)) * prod(resolution);
named_imgs_volumes(EC_translated) = sum(img(:)) * voxel_volume;
all_compartment_volumes(1) = named_imgs_volumes(EC_translated);

all_compartment_indices = containers.Map('KeyType', 'char', 'ValueType', 'double');
for i = 1:n_compartments
    all_compartment_indices(names_sorted{i}) = i;
end
%}


% Create indexed image

% Combine first
[all_compartments_image, named_imgs, names_sorted] = combine_compartment_images(named_imgs, names_sorted, voxel_volume, network_info.compartments);

all_compartments_image_before_downsampling = all_compartments_image;
named_imgs_before_downsampling = named_imgs;
resolution_before_downsampling = resolution;
resolution_single_before_downsampling = resolution_single;
voxel_volume_before_downsampling = voxel_volume;


if any(downsampling ~= 1)
    combine_resize_compartment_images_result = combine_resize_compartment_images(named_imgs, voxel_volume, names_sorted, network_info.compartments, ones(1, 3) .* downsampling);
    all_compartments_image = combine_resize_compartment_images_result.all_compartments_image;
    named_imgs = combine_resize_compartment_images_result.named_imgs;
    actual_scale = combine_resize_compartment_images_result.actual_scale;
    named_imgs_exclusive_volumes = combine_resize_compartment_images_result.named_imgs_exclusive_volumes;
    names_sorted = combine_resize_compartment_images_result.names_sorted;

    resolution = resolution ./ actual_scale;
    resolution_single = resolution_single ./ power(prod(actual_scale), 1/3);
    voxel_volume = voxel_volume ./ prod(actual_scale);

    for i = 1:length(names)
        name = names{i};
        img = named_imgs(name);
        % named_imgs_volumes(name) = sum(img(:)) * prod(resolution);
        named_imgs_volumes(name) = sum(img(:)) * voxel_volume;
    end
else
    [named_imgs, named_imgs_exclusive_volumes] = get_exclusive_images(all_compartments_image, voxel_volume, names_sorted);
end


all_compartment_volumes = zeros(n_compartments, 1);
for i = 1:n_compartments
    all_compartment_volumes(i) = named_imgs_volumes(names_sorted{i});
end

% Exclusive volumes (do not include volumes of compartments inside)
all_compartment_exclusive_volumes = zeros(n_compartments, 1);
for i = 1:n_compartments
    all_compartment_exclusive_volumes(i) = named_imgs_exclusive_volumes(names_sorted{i});
end

if options.output.indexedimage
    savepath = options.output.temp.savepath;
    if any(downsampling ~= 1)
        imwrite(reshape_contrast(single(all_compartments_image_before_downsampling), -1), [savepath, ' all_compartments_image_before_downsampling.png']);
    end
    imwrite(reshape_contrast(single(all_compartments_image), -1), [savepath, ' all_compartments_image.png']);
end


% Assign volumes and surface areas in network_info.compartments
for i = 1:length(named_imgs)
    name = names_sorted{i};
    compartment = network_info.compartments(name);
    compartment.size_expression = named_imgs_exclusive_volumes(name);
    network_info.compartments(name) = compartment;
end


% Find adjacent compartments before downsampling
if use_image_adjacency
    % Find adjacent compartments
    adjacent_pairs = getAdjacentValues(all_compartments_image, connectivity);
else
    adjacent_pairs = getAdjacentCompartments(network_info.compartments, all_compartment_indices);
    
    adjacent_pairs = sort(adjacent_pairs, 2);
    adjacent_pairs = unique(adjacent_pairs, 'rows');
end


% Print for use in BNGL files
fprintf('\n');
fprintf('\n');
fprintf('Compartment volumes (exclusive):\n');
for i = 1:n_compartments
    fprintf('    %s: %.6f * %s\n', names_sorted{i}, named_imgs_exclusive_volumes(names_sorted{i}), char(output_volume_unit));
end
fprintf('\n');
fprintf('Compartment volumes (inclusive):\n');
for i = 1:n_compartments
    fprintf('    %s: %.6f * %s\n', names_sorted{i}, named_imgs_volumes(names_sorted{i}), char(output_volume_unit));
end
fprintf('\n');
fprintf('\n');


% Add indices to network_info.compartments
n_network_info_compartments = length(network_info.compartments);
network_info_compartments_keys = network_info.compartments.keys;
network_info_compartments_keys_indices_map = containers.Map('KeyType', 'char', 'ValueType', 'double');
for i = 1:n_network_info_compartments
    compartment_name = network_info_compartments_keys{i};
    network_info_compartments_keys_indices_map(compartment_name) = i;
    compartment = network_info.compartments(compartment_name);
    compartment.image_index = nan;
    compartment.image_outside_index = nan;
    compartment.image_insides_indices = [];
    compartment.insides = {};
    network_info.compartments(compartment_name) = compartment;
end

for i = 1:n_network_info_compartments
    compartment_name = network_info_compartments_keys{i};
    compartment = network_info.compartments(compartment_name);
    compartment_outside_name = compartment.outside;
    if all_compartment_indices.isKey(compartment_name)
        compartment.image_index = all_compartment_indices(compartment_name);
    end
    if all_compartment_indices.isKey(compartment_outside_name)
        compartment.image_outside_index = all_compartment_indices(compartment_outside_name);
        compartment_outside = network_info.compartments(compartment_outside_name);
        compartment_outside.image_insides_indices(end+1) = i;
        network_info.compartments(compartment_outside_name) = compartment_outside;
    end
    if ~isempty(compartment_outside_name)
        compartment_outside = network_info.compartments(compartment_outside_name);
        compartment_outside.insides{end+1} = compartment_name;
        network_info.compartments(compartment_outside_name) = compartment_outside;
    end
    network_info.compartments(compartment_name) = compartment;
end

% Create compartment adjacency matrices for paths of length 1 to adjacent_pairs_max_degree
network_info_adjacency_max_degree = 2;
network_info_adjacency_matrix = false(n_network_info_compartments, n_network_info_compartments);
for i = 1:n_network_info_compartments
    compartment_name = network_info_compartments_keys{i};
    compartment = network_info.compartments(compartment_name);
    compartment_insides = compartment.insides;
    compartment_insides_indices = cellfun(@(x)network_info_compartments_keys_indices_map(x), compartment.insides);
    network_info_adjacency_matrix(i, compartment_insides_indices) = true;
    network_info_adjacency_matrix(compartment_insides_indices, i) = true;
end

network_info_adjacency_matrices = repmat(network_info_adjacency_matrix, 1, 1, network_info_adjacency_max_degree);
for i = 2:network_info_adjacency_max_degree
    network_info_adjacency_matrices(:, :, i) = single(network_info_adjacency_matrices(:, :, i-1)) * network_info_adjacency_matrix > 0;
end

function result = is_network_info_compartment_pair_adjacent(given_compartment1, given_compartment2, given_degrees)
    if nargin < 3
        given_degrees = 1;
    end
    given_compartment1_index = given_compartment1;
    given_compartment2_index = given_compartment2;
    if ischar(given_compartment1_index)
        given_compartment1_index = network_info_compartments_keys_indices_map(given_compartment1_index);
    end
    if ischar(given_compartment2_index)
        given_compartment2_index = network_info_compartments_keys_indices_map(given_compartment2_index);
    end
    result = network_info_adjacency_matrices(given_compartment1_index, given_compartment2_index, given_degrees);
end


% Create compartment adjacency matrices for paths of length 1 to adjacent_pairs_max_degree
adjacent_pairs_max_degree = 2;
adjacent_pairs_matrix = false(n_compartments, n_compartments);
adjacent_pairs_matrices = repmat(adjacent_pairs_matrix, 1, 1, adjacent_pairs_max_degree);
for i = 2:adjacent_pairs_max_degree
    adjacent_pairs_matrices(:, :, i) = single(adjacent_pairs_matrices(:, :, i-1)) * adjacent_pairs_matrix > 0;
end
function result = is_compartment_pair_adjacent(given_compartment1, given_compartment2, given_degrees)
    if nargin < 3
        given_degrees = 1;
    end
    given_compartment1_index = given_compartment1;
    given_compartment2_index = given_compartment2;
    if ischar(given_compartment1_index)
        given_compartment1_index = all_compartment_indices(given_compartment1_index);
    end
    if ischar(given_compartment2_index)
        given_compartment2_index = all_compartment_indices(given_compartment2_index);
    end
    result = adjacent_pairs_matrices(given_compartment1_index, given_compartment2_index, given_degrees);
end


% Find all objects in each compartment

all_compartments_objects = containers.Map('KeyType', 'char', 'ValueType', 'any');
% New indexed image for individual compartmental objects
% all_compartments_objects_image_before_downsampling = [];
all_compartments_objects_image = [];
all_compartments_objects_names = {};
all_compartments_objects_names_to_compartment_names = containers.Map('KeyType', 'char', 'ValueType', 'char');
all_compartments_objects_indices = containers.Map('KeyType', 'char', 'ValueType', 'double');
should_only_keep_largest_framework_objects = true;

% use_image_adjacency_for_compartments = false;
use_image_adjacency_for_compartments = true;

if use_image_adjacency_for_compartments
    all_compartments_objects_start_time = tic();

    xyz_to_ijk = [2, 1, 3];
    ijk_to_xyz = xyz_to_ijk;

    k = 1;
    for i = 1:n_compartments
        name = names_sorted{i};
        name_img = named_imgs(name);
        name_img_exclusive = all_compartments_image == i;
        name_img_objects = bwconncomp(name_img_exclusive, connectivity);
        if strcmp(name, EC_translated)
            if i ~= 1
                error('"EC"/"%s" should be first in names_sorted', EC_translated);
            end
            % EC can be split into multiple objects, assume it should be considered one
            name_img_objects.PixelIdxList = {cell2mat(name_img_objects.PixelIdxList')};
            name_img_objects.NumObjects = 1;
        end
        
        name_img_objects_stats = regionprops(name_img_objects);
        
        should_only_keep_largest_compartment_objects = should_only_keep_largest_framework_objects && is_framework_compartment_name_function(name) && length(name_img_objects_stats) > 1;
        
        if should_only_keep_largest_compartment_objects
            compartment_largest_object_index = nan;
            compartment_largest_object_volume = 0;
        end
        
        compartment_mask = all_compartments_image == i;
        previous_all_compartments_objects_image = all_compartments_objects_image;
        for j = 1:length(name_img_objects_stats)
            compartment_object_name = [name, num2str(k)];
            temp = name_img_objects_stats(j);
            [temp_pixels_i, temp_pixels_j, temp_pixels_k] = ind2sub(size(all_compartments_image), name_img_objects.PixelIdxList{j});
            temp_pixels_ijk = [temp_pixels_i, temp_pixels_j, temp_pixels_k];
            temp_pixels_xyz = temp_pixels_ijk(:, ijk_to_xyz);
            temp_bounding_box_2d = reshape(temp.BoundingBox, [], 2)';
            temp_lower_corner_xyz = ceil(temp.BoundingBox(1:3));
            temp_lower_corner_ijk = temp_lower_corner_xyz(xyz_to_ijk);
            temp_size_xyz = temp.BoundingBox(4:6);
            temp_size_ijk = temp_size_xyz(xyz_to_ijk);
            temp_upper_corner_xyz = temp_lower_corner_xyz + temp_size_xyz - 1;
            temp_upper_corner_ijk = temp_upper_corner_xyz(xyz_to_ijk);
            temp_crop_indices_cell = {temp_lower_corner_ijk(1):temp_upper_corner_ijk(1), temp_lower_corner_ijk(2):temp_upper_corner_ijk(2), temp_lower_corner_ijk(3):temp_upper_corner_ijk(3)};
            temp_pixels_xyz_cropped = temp_pixels_xyz - repmat(temp_lower_corner_xyz - 1, size(temp_pixels_xyz, 1), 1);
            temp_pixels_ijk_cropped = temp_pixels_xyz_cropped(:, xyz_to_ijk);
            temp_pixels_ind_cropped = sub2ind(temp_size_ijk, temp_pixels_ijk_cropped(:, 1), temp_pixels_ijk_cropped(:, 2), temp_pixels_ijk_cropped(:, 3));
            mask = false(temp_size_ijk);
            mask(temp_pixels_ind_cropped) = true;
            % mask = mask & compartment_mask(temp_crop_indices_cell{:});
            
            temp.SamplePixel = temp_pixels_xyz(1, :);
            mask_filled = imfill(mask, connectivity, 'holes');
            
            mask_eroded = imerode(mask, connectivity_se);
            mask_filled_eroded = imerode(mask_filled, connectivity_se);
            
            mask_sum = sum(mask(:));
            mask_filled_sum = sum(mask_filled(:));
            mask_eroded_sum = sum(mask_eroded(:));
            mask_filled_eroded_sum = sum(mask_filled_eroded(:));
            
            temp.Volume = mask_sum;
            temp.FilledVolume = mask_filled_sum;
            temp.Area = mask_sum - mask_eroded_sum;
            temp.FilledArea = mask_filled_sum - mask_filled_eroded_sum;
            
            % TODO: Correctly compute area of scaled geometry instead of approximating it by scaling voxel counts likes this
            temp.Volume = temp.Volume * voxel_volume;
            temp.FilledVolume = temp.FilledVolume * voxel_volume;
            temp.Area = temp.Area * voxel_face_area;
            temp.FilledArea = temp.FilledArea * voxel_face_area;
            
            % if temp.FilledVolume == numel(mask)
            if strcmp(name, EC_translated)
                temp.FilledVolumeUnlessAll = temp.Volume;
                temp.FilledAreaUnlessAll = temp.Area;
            else
                temp.FilledVolumeUnlessAll = temp.FilledVolume;
                temp.FilledAreaUnlessAll = temp.FilledArea;
            end
            
            temp.name = compartment_object_name;
            temp.compartment_name = name;
            temp.index = k;
            temp.subVolumeName = name;
            temp.SpatialObjectName = ['vobj_', compartment_object_name];
        
            if should_only_keep_largest_compartment_objects && compartment_largest_object_volume < temp.Volume
                if ~isnan(compartment_largest_object_index)
                    % Remove previous largest object
                    compartment_object_to_remove_name = all_compartments_objects_names{compartment_largest_object_index};
                    compartment_object_to_remove = all_compartments_objects(compartment_object_to_remove_name);
                    all_compartments_objects.remove(compartment_object_to_remove_name);
                    all_compartments_objects_names_to_compartment_names.remove(compartment_object_to_remove_name);
                    all_compartments_objects_indices.remove(compartment_object_to_remove_name);
                    old_mask = all_compartments_objects_image == compartment_largest_object_index;
                    all_compartments_objects_image(old_mask) = previous_all_compartments_objects_image(old_mask);
                    all_compartments_objects_names = all_compartments_objects_names(1:end - 1);
                    k = k - 1;
                end
                compartment_largest_object_index = temp.index;
                compartment_largest_object_volume = temp.Volume;
            end
        
            if ~should_only_keep_largest_compartment_objects || compartment_largest_object_index == temp.index
                all_compartments_objects(compartment_object_name) = temp;
                
                % Image
                if isempty(all_compartments_objects_image)
                    all_compartments_objects_image = zeros(size(all_compartments_image), 'uint32');
                end
                
                all_compartments_objects_image_cropped = all_compartments_objects_image(temp_crop_indices_cell{:});
                all_compartments_objects_image_cropped(mask) = k;
                all_compartments_objects_image(temp_crop_indices_cell{:}) = all_compartments_objects_image_cropped;
                
                all_compartments_objects_names{k} = compartment_object_name;
                all_compartments_objects_names_to_compartment_names(compartment_object_name) = name;
                all_compartments_objects_indices(compartment_object_name) = k;
                
                k = k + 1;
            end
        end
    end

    n_compartments_objects = length(all_compartments_objects_names);

    all_compartments_objects_elapsed_time = toc(all_compartments_objects_start_time);

else

    % Generate meshes from CSGdata and add to meshData
    
    if ~should_only_keep_largest_framework_objects
        error('Not implemented');
    end

    all_compartments_objects_start_time = tic();

    % options.output.name_map = name_map;
    % options.output.resolution_before_downsampling = resolution_before_downsampling;
    options.output.resolution_single_before_downsampling = resolution_single_before_downsampling;
    % options.output.model_names = model_names;
    
    % [CSGdata, meshData] = convertCSGToMesh(CSGdata, meshData, models, imgs, options);
    
    %{
    % Rename meshes named for membranes after the volumes they contain
    for i = 1:length(meshData)
        meshData(i).source = 'mesh';
        for j = 1:length(meshData(i).list)
            meshData_list_item_name = meshData(i).list(j).name;
            meshData_list_item_name = translateWithDefaultIdentity(name_map, meshData_list_item_name);
            % if ~is_membrane_function(meshData_list_item_name)
            if network_info.compartments.isKey(meshData_list_item_name) && network_info.compartments(meshData_list_item_name).spatial_dimensions == 2
                % meshData_list_item_name = network_info.compartments(meshData_list_item_name).outside;
                meshData_list_item_name = network_info.compartments(meshData_list_item_name).insides{1};
            end
            meshData(i).list(j).name = meshData_list_item_name;
        end
    end
    %}
    
    % Flatten objects into all_compartments_objects and add geometric and relational data
    
    warning('CellOrganizer:readNetworkIntoGeometry', 'No implemented difference between mesh Volume and FilledVolume, etc.');
    k = 0;
    framework_comparment_object_names = containers.Map('KeyType', 'char', 'ValueType', 'char');
    for i = 1:length(meshData)
        single_meshData = meshData(i);
        single_meshData_name = single_meshData.name;
        single_meshData_list = single_meshData.list;
        for j = 1:length(single_meshData_list)
            k = k + 1;
            single_meshData_list_item = single_meshData_list(j);
            single_meshData_list_item_name = single_meshData_list_item.name;
            single_meshData_list_item_mesh = single_meshData_list_item.mesh;
            compartment_name = single_meshData_list_item_name;
            compartment_index = all_compartment_indices(compartment_name);
            compartment_object_name = [single_meshData_list_item_name, num2str(k-1)];
            if is_framework_compartment_name_function(compartment_name)
                framework_comparment_object_names(compartment_name) = compartment_object_name;
            end
            
            
            temp = struct();
            temp.BoundingBox = [min(single_meshData_list_item_mesh.vertices, [], 1), max(single_meshData_list_item_mesh.vertices, [], 1)];
            
            single_meshData_list_item_mesh = single_meshData_list_item_mesh;
            single_meshData_list_item_mesh_statistics = triangleMeshStatistics(single_meshData_list_item_mesh);
            
            temp.Volume = nan;
            temp.FilledVolume = single_meshData_list_item_mesh_statistics.volume;
            temp.Area = nan;
            temp.FilledArea = single_meshData_list_item_mesh_statistics.area;
            
            if strcmp(single_meshData_list_item_name, EC_translated)
                temp.FilledVolumeUnlessAll = temp.Volume;
                temp.FilledAreaUnlessAll = temp.Area;
            else
                temp.FilledVolumeUnlessAll = temp.FilledVolume;
                temp.FilledAreaUnlessAll = temp.FilledArea;
            end
            
            temp.name = compartment_object_name;
            temp.compartment_name = compartment_name;
            temp.subVolumeName = compartment_name;
            temp.SpatialObjectName = ['vobj_', compartment_object_name];
            temp.index = k;
            
            all_compartments_objects(compartment_object_name) = temp;
            all_compartments_objects_names{end+1} = temp.name;
            all_compartments_objects_names_to_compartment_names(compartment_object_name) = temp.compartment_name;
            all_compartments_objects_indices(compartment_object_name) = temp.index;
        end
    end
    
    n_compartments_objects = length(all_compartments_objects_names);
    
    % Subtract interior volumes
    for i = 1:n_compartments_objects
        object_name = all_compartments_objects_names{i};
        object = all_compartments_objects(object_name);
        object.Volume = object.FilledVolume;
        object.Area = object.FilledArea;
        all_compartments_objects(object_name) = object;
    end
    for i = 1:n_compartments_objects
        object_name = all_compartments_objects_names{i};
        object = all_compartments_objects(object_name);
        object_compartment_name = object.compartment_name;
        % if compartment_object.FilledVolume == numel(mask)
        if strcmp(object_compartment_name, EC_translated)
            object.FilledVolumeUnlessAll = object.Volume;
            object.FilledAreaUnlessAll = object.Area;
        else
            object.FilledVolumeUnlessAll = object.FilledVolume;
            object.FilledAreaUnlessAll = object.FilledArea;
        end
        all_compartments_objects(object_name) = object;
        % Modify outside compartment object's volume and area
        object_compartment = network_info.compartments(object_compartment_name);
        object_outside_compartment_name = object_compartment.double_outside;
        if ~isempty(object_outside_compartment_name)
            object_outside_compartment_object_name = framework_comparment_object_names(object_outside_compartment_name);
            object_outside_compartment_object = all_compartments_objects(object_outside_compartment_object_name);
            object_outside_compartment_object.Volume = object_outside_compartment_object.Volume - object.FilledVolume;
            object_outside_compartment_object.Area = object_outside_compartment_object.Area + object.FilledArea;
            all_compartments_objects(object_outside_compartment_object_name) = object_outside_compartment_object;
        end
    end

end



% Assign insides and outsides for objects

compartments_compartment_objects_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
all_compartments_objects_keys = all_compartments_objects.keys();
for i = 1:length(all_compartments_objects_keys)
    compartment_object_name = all_compartments_objects_keys{i};
    compartment_object = all_compartments_objects(compartment_object_name);
    compartment_name = compartment_object.compartment_name;
    if ~compartments_compartment_objects_map.isKey(compartment_name)
        compartments_compartment_objects_map(compartment_name) = {};
    end
    compartments_compartment_objects_map(compartment_name) = [compartments_compartment_objects_map(compartment_name), {compartment_object_name}];
    n_compartment_objects = length(compartments_compartment_objects_map(compartment_name));
end

empty_compartment = struct('double_outside', '', 'double_insides', {{}});
for i = 1:length(all_compartments_objects_keys)
    compartment_object_name = all_compartments_objects_keys{i};
    compartment_object = all_compartments_objects(compartment_object_name);
    compartment_name = compartment_object.compartment_name;
    if network_info.compartments.isKey(compartment_name)
        compartment = network_info.compartments(compartment_name);
    else
        compartment = empty_compartment;
    end
    % compartment_outside = compartment.double_outside;
    % compartment_insides = compartment.double_insides;
    compartment_outside = compartment.double_outside;
    compartment_insides = compartment.double_insides;
    compartment_object_compartment_object_outside = '';
    if compartments_compartment_objects_map.isKey(compartment_outside)
        compartment_object_compartment_object_outside = compartments_compartment_objects_map(compartment_outside);
        % This assumes surface compartments can be contained only by the first of other surface compartments (for example, geometry for only one cell containing multiple vesicles is supported)
        compartment_object_compartment_object_outside = compartment_object_compartment_object_outside{1};
    end
    compartment_object_compartment_objects_insides = {};
    for j = 1:length(compartment_insides)
        if compartments_compartment_objects_map.isKey(compartment_insides{j})
            compartment_object_compartment_objects_insides = [compartment_object_compartment_objects_insides, compartments_compartment_objects_map(compartment_insides{j})];
        end
    end
    compartment_object.insides = compartment_object_compartment_objects_insides;
    compartment_object.outside = compartment_object_compartment_object_outside;
    all_compartments_objects(compartment_object_name) = compartment_object;
end


n_compartments_objects = length(all_compartments_objects_names);



if use_image_adjacency
    % Find adjacent compartment objects
    object_adjacent_pairs_start_time = tic();
    object_adjacent_pairs = getAdjacentValues(all_compartments_objects_image, connectivity);
    % object_adjacent_pairs_before_downsampling = object_adjacent_pairs;
    object_adjacent_pairs_elapsed_time = toc(object_adjacent_pairs_start_time);

    warning('CellOrganizer:readNetworkIntoGeometry', 'Some inconsistencies! Voxel not 6-connected then determined to be 6-adjacent')
else
    object_adjacent_pairs = getAdjacentObjects(all_compartments_objects);
    % object_adjacent_pairs_before_downsampling = object_adjacent_pairs;
end




% Following based on a VCML file saved using Virtual Cell



% Data required to generate VCML

all_compartment_data = struct();
for i = 1:n_compartments
    name = names_sorted{i};
    compartment = network_info.compartments(name);
    object_volume = all_compartment_exclusive_volumes(i);
    all_compartment_data(i).name = name;
    all_compartment_data(i).object_volume = object_volume;
    all_compartment_data(i).FeatureName = name;
    all_compartment_data(i).DiagramNodeName = all_compartment_data(i).FeatureName;
    all_compartment_data(i).DiagramName = all_compartment_data(i).DiagramNodeName;
    all_compartment_data(i).PixelClassNodeName = name;
    all_compartment_data(i).SubVolumeNodeName = name;
    all_compartment_data(i).double_outside = compartment.double_outside;
    all_compartment_data(i).double_insides = compartment.double_insides;
end


all_object_data = struct();
for i = 1:n_compartments_objects
    name = all_compartments_objects_names{i};
    object = all_compartments_objects(name);
    object_volume = object.Volume;
    all_object_data(i).name = name;
    all_object_data(i).object = object;
    all_object_data(i).object_volume = object_volume;
    all_object_data(i).SpatialObjectName = object.SpatialObjectName;
    all_object_data(i).VolumeRegionNodeName = name;
    all_object_data(i).SubVolumeNodeName = all_compartments_objects_names_to_compartment_names(name);
end

all_object_membrane_data = struct();
k = 0;
for ij = object_adjacent_pairs'
    k = k + 1;
    i = ij(1);
    j = ij(2);
    name_i = all_compartments_objects_names{i};
    name_j = all_compartments_objects_names{j};
    name = [name_i, '_', name_j];
    object_i = all_compartments_objects(name_i);
    object_j = all_compartments_objects(name_j);
    object_i_volume = object_i.FilledVolumeUnlessAll;
    object_j_volume = object_j.FilledVolumeUnlessAll;
    object_i_area = object_i.FilledAreaUnlessAll;
    object_j_area = object_j.FilledAreaUnlessAll;
    % Assumes one compartment surrounds the other
    membrane_area = min(object_i_area, object_j_area);
    all_object_membrane_data(k).ij = ij;
    all_object_membrane_data(k).i = i;
    all_object_membrane_data(k).j = j;
    all_object_membrane_data(k).name = name;
    all_object_membrane_data(k).name_i = name_i;
    all_object_membrane_data(k).name_j = name_j;
    all_object_membrane_data(k).object_i = object_i;
    all_object_membrane_data(k).object_j = object_j;
    all_object_membrane_data(k).object_i_volume = object_i_volume;
    all_object_membrane_data(k).object_j_volume = object_j_volume;
    all_object_membrane_data(k).object_i_area = object_i_area;
    all_object_membrane_data(k).object_j_area = object_j_area;
    all_object_membrane_data(k).object_area = membrane_area;
    all_object_membrane_data(k).MembraneRegionNodeName = ['membrane_', name];
    all_object_membrane_data(k).SpatialObjectName = ['sobj_', name];
    if object_i_volume > object_j_volume
        all_object_membrane_data(k).subVolumeOutside = all_compartments_objects_names_to_compartment_names(name_i);
        all_object_membrane_data(k).subVolumeInside = all_compartments_objects_names_to_compartment_names(name_j);
        all_object_membrane_data(k).regionIdOutside = i-1;
        all_object_membrane_data(k).regionIdInside = j-1;
    else
        all_object_membrane_data(k).subVolumeOutside = all_compartments_objects_names_to_compartment_names(name_j);
        all_object_membrane_data(k).subVolumeInside = all_compartments_objects_names_to_compartment_names(name_i);
        all_object_membrane_data(k).regionIdOutside = j-1;
        all_object_membrane_data(k).regionIdInside = i-1;
    end
end


all_membrane_data_start_time = tic();

all_membrane_data = struct('i', {}, 'j', {}, 'outside_object', {}, 'inside_object', {}, 'pair_string', {}, 'name', {}, 'area', {}, 'MembraneNodeName', {}, 'MembraneVoltageName', {}, 'DiagramNodeName', {}, 'SurfaceClassNodeName', {}, 'SubVolume1Ref', {}, 'SubVolume2Ref', {}, 'OutsideCompartment', {}, 'InsideCompartment', {});

k = 1;
ball_se = connectivity_se;
for ij = adjacent_pairs'
%{
% Ignore topological changes produced by downsampling
for ij = adjacent_pairs_before_downsampling'
%}
    i = ij(1);
    j = ij(2);
    start_obj = all_compartment_data(i);
    link_obj = all_compartment_data(j);
    
    if strcmp(start_obj.name, link_obj.double_outside)
        outside_obj = start_obj;
        inside_obj = link_obj;
    else
        outside_obj = link_obj;
        inside_obj = start_obj;
    end
    
    % Sort to make independent of order
    [pair_string_translated, pair_string_cell, pair_string] = get_compartment_pair_string(start_obj.name, link_obj.name, name_map);
    name = pair_string_translated;
    name_membrane = [name, '_membrane'];
    
    % Estimate surface area from dilation overlap
    start_image_dilated = single(imdilate(all_compartments_image == i, ball_se));
    link_image_dilated = single(imdilate(all_compartments_image == j, ball_se));
    area = start_image_dilated & link_image_dilated;
    area = sum(area(:)) / 2;
    area = area * voxel_face_area;
    
    % Compute surface area from boundaries. Assumes cubic voxels.
    temp_img_i = all_compartments_image == i;
    temp_img_j = all_compartments_image == j;
    area = 0;
    % all_compartments_image changing from i to j in positive y, x, and z directions
    area = area + sum(reshape(temp_img_i(1:end-1, :, :) & temp_img_j(2:end, :, :), 1, []));
    area = area + sum(reshape(temp_img_i(:, 1:end-1, :) & temp_img_j(:, 2:end, :), 1, []));
    area = area + sum(reshape(temp_img_i(:, :, 1:end-1) & temp_img_j(:, :, 2:end), 1, []));
    % all_compartments_image changing from j to i in positive y, x, and z directions
    area = area + sum(reshape(temp_img_j(1:end-1, :, :) & temp_img_i(2:end, :, :), 1, []));
    area = area + sum(reshape(temp_img_j(:, 1:end-1, :) & temp_img_i(:, 2:end, :), 1, []));
    area = area + sum(reshape(temp_img_j(:, :, 1:end-1) & temp_img_i(:, :, 2:end), 1, []));
    area = area * voxel_face_area;
    
    all_membrane_data(k).i = i;
    all_membrane_data(k).j = j;
    all_membrane_data(k).outside_object = outside_obj;
    all_membrane_data(k).inside_object = inside_obj;
    all_membrane_data(k).pair_string = pair_string;
    all_membrane_data(k).name = name;
    all_membrane_data(k).area = area;
    all_membrane_data(k).MembraneNodeName = name;
    all_membrane_data(k).MembraneVoltageName = ['Voltage_', name];
    all_membrane_data(k).DiagramNodeName = name;
    all_membrane_data(k).MembraneSubDomainNodeName = name_membrane;
    all_membrane_data(k).SurfaceClassNodeName = name_membrane;
    all_membrane_data(k).SubVolume1Ref = start_obj.name;
    all_membrane_data(k).SubVolume2Ref = link_obj.name;
    all_membrane_data(k).OutsideCompartment = outside_obj.name;
    all_membrane_data(k).InsideCompartment = inside_obj.name;
    
    k = k + 1;
end
all_membrane_names = {all_membrane_data.name};
function result = is_membrane_function(given_compartments)
    if ischar(given_compartments)
        given_compartments = {given_compartments};
    end
    result = cellfun(@(x)any(strcmp(x, all_membrane_names)), given_compartments);
end

all_membrane_data_elapsed_time = toc(all_membrane_data_start_time);


% Print for use in BNGL files
fprintf('\n');
fprintf('\n');
fprintf('Membrane areas:\n');
for k = 1:length(all_membrane_data)
    object = all_membrane_data(k);
    name = object.name;
    area = object.area;
    fprintf('    %s: %.6f * %s\n', name, area, char(output_area_unit));
end
fprintf('\n');
fprintf('\n');






% Add species names
species_digits = 0;
species_name_format = '';
function result = species_index_to_name_function(given_species_index)
    result = sprintf(species_name_format, given_species_index);
end
function updateSpeciesNames()
    maximum_number_species = length(network_info.species);
    if add_translocation_intermediates
        % Include maximum number of intermediate species
        maximum_number_species = maximum_number_species * length(network_info.compartments);
    end
    species_digits = ceil(log10(maximum_number_species+1));
    species_name_format = ['s%0', num2str(species_digits), 'i'];
    species_names = cellfun(@species_index_to_name_function, num2cell(1:length(network_info.species)), 'UniformOutput', false);
    species_extended_names = cellfun(@(x)network_info.species(x).species_graph, num2cell(1:length(network_info.species)), 'UniformOutput', false);
    [network_info.species.name] = species_names{:};
    [network_info.species.extended_name] = species_extended_names{:};
end
updateSpeciesNames();



% Set or infer parameter units
warning('CellOrganizer:readNetworkIntoGeometry', 'Assuming all unknown parameters have no units');
parameters_names_backup_map = containers.Map('KeyType', 'char', 'ValueType', 'char');
for i = 1:length(network_info_parameters_keys)
    parameter_name = network_info_parameters_keys{i};
    parameter = network_info.parameters(parameter_name);
    % parameter_value_expression = parameter.value_expression;
    parameter_units = [];
    if strbegins(parameter_name, compartment_volume_prefix) || strends(parameter_name, '_volume')
        parameter_units = net_volume_unit;
    elseif strbegins(parameter_name, compartment_area_prefix) || strends(parameter_name, '_area')
        parameter_units = net_area_unit;
    elseif strbegins(parameter_name, compartment_diffusion_coefficient_prefix)
        parameter_units = net_diffusion_coefficient_unit;
    elseif strends(parameter_name, '_length') || strcmp(parameter_name, 'eff_width')
        parameter_units = net_length_unit;
    elseif ~isempty(regexp(parameter_name, '^_rateLaw[0-9]+$'))
        warning('CellOrganizer:readNetworkIntoGeometry', 'Not setting _rateLaw* parameter dimensions');
    elseif strends(parameter_name, '_conc')
        parameter_units = net_concentration_unit;
    elseif ~isempty(regexp(parameter_name, '^_InitialConc[0-9]+$')) || strends(parameter_name, '_count')
        parameter_units = net_count_unit;
    end
    if ~isempty(parameter_units)
        original_parameter = parameter;
        parameter = DimensionedExpression(parameter.value, parameter_units);
        updateParametersMap(parameter_name, parameter);
        parameter_to_add = DimensionedExpression(parameter);
        parameter_to_add_name = [backup_prefix, parameter_name];
        updateParametersMap(parameter_to_add_name, parameter_to_add);
        parameters_names_backup_map(parameter_name) = parameter_to_add_name;
    end
end
% Replace known values
updateParametersMap('deci', DimensionedExpression(1e-1));
updateParametersMap('centi', DimensionedExpression(1e-2));
updateParametersMap('milli', DimensionedExpression(1e-3));
updateParametersMap('micro', DimensionedExpression(1e-6));
updateParametersMap('nano', DimensionedExpression(1e-9));
updateParametersMap('pico', DimensionedExpression(1e-12));
updateParametersMap('femto', DimensionedExpression(1e-15));
updateParametersMap('second', DimensionedExpression(1, 's'));
updateParametersMap('molecule', DimensionedExpression(1, 'molecule'));
updateParametersMap('meter', DimensionedExpression(1, 'm'));
updateParametersMap('meter2', parameterReference('meter')^2);
updateParametersMap('meter3', parameterReference('meter')^3);
updateParametersMap('centimeter', DimensionedExpression(1, 'cm'));
updateParametersMap('centimeter2', parameterReference('centimeter')^2);
updateParametersMap('centimeter3', parameterReference('centimeter')^3);
updateParametersMap('micrometer', DimensionedExpression(1, 'um'));
updateParametersMap('micrometer2', parameterReference('micrometer')^2);
updateParametersMap('micrometer3', parameterReference('micrometer')^3);
updateParametersMap('micron', parameterReference('micrometer')^1);
updateParametersMap('micron2', parameterReference('micron')^2);
updateParametersMap('micron3', parameterReference('micron')^3);
updateParametersMap('nanometer', DimensionedExpression(1, 'nm'));
updateParametersMap('nanometer2', parameterReference('nanometer')^2);
updateParametersMap('nanometer3', parameterReference('nanometer')^3);
updateParametersMap('angstrom', DimensionedExpression(1e2, 'pm'));
updateParametersMap('angstrom2', parameterReference('angstrom')^2);
updateParametersMap('angstrom3', parameterReference('angstrom')^3);
updateParametersMap('picometer', DimensionedExpression(1, 'nm'));
updateParametersMap('picometer2', parameterReference('picometer')^2);
updateParametersMap('picometer3', parameterReference('picometer')^3);
updateParametersMap('mole', DimensionedExpression(1, 'mol'));
updateParametersMap('N_Avo', DimensionedExpression('6.022140857e23', 'molecule.mol-1'));
updateParametersMap('liter', DimensionedExpression(1, 'L'));
updateParametersMap('mole', DimensionedExpression(1, 'mol'));
updateParametersMap('molar', DimensionedExpression(1, 'M'));
updateParametersMap('millimolar', DimensionedExpression(1, 'mM'));
updateParametersMap('micromolar', DimensionedExpression(1, 'uM'));
updateParametersMap('nanomolar', DimensionedExpression(1, 'nM'));
updateParametersMap('picomolar', DimensionedExpression(1, 'pM'));
updateParametersMap('concentration', DimensionedExpression(1, net_concentration_unit));

% updateParametersMap();
inferUnits();



% Set diffusion coefficient units

warning('CellOrganizer:readNetworkIntoGeometry', 'Assuming dc_* parameters are diffusion coefficients');
network_info_species_diffusion_coefficients = repmat({nan}, size(network_info.species));
[network_info.species(:).diffusion_coefficient] = network_info_species_diffusion_coefficients{:};
for i = 1:length(network_info.species)
    single_species = network_info.species(i);
    name = single_species.name;
    single_species_compartments = single_species.compartments;
    diffusion_coefficient_parameter_name = [compartment_diffusion_coefficient_prefix, name];
    if ~network_info.parameters.isKey(diffusion_coefficient_parameter_name)
        if isnan(single_species.diffusion_coefficient)
            % Default overall
            single_species_compartments_diffusion_coefficients = repmat(output_diffusion_coefficient, size(single_species_compartments));
            % Default for compartments
            for j = 1:length(single_species_compartments)
                compartment_name = single_species_compartments{j};
                compartment_diffusion_coefficient_parameter_name = [compartment_diffusion_coefficient_prefix, compartment_name];
                if network_info.parameters.isKey(compartment_diffusion_coefficient_parameter_name)
                    single_species_compartments_diffusion_coefficients(j) = parameterReference(compartment_diffusion_coefficient_parameter_name);
                end
            end
            % Geometric mean
            single_species_diffusion_coefficient = prod(single_species_compartments_diffusion_coefficients);
            single_species_diffusion_coefficient = single_species_diffusion_coefficient^(1 / length(single_species_compartments_diffusion_coefficients));
            single_species.diffusion_coefficient = single_species_diffusion_coefficient;
        end
        % network_info.parameters(diffusion_coefficient_parameter_name) = single_species.diffusion_coefficient;
        updateParametersMap(diffusion_coefficient_parameter_name, single_species.diffusion_coefficient);
    end
    single_species.diffusion_coefficient = parameterReference(diffusion_coefficient_parameter_name);
    network_info.species(i) = single_species;
end

% updateParametersMap();
inferUnits();



% Replace references to backed up volume and surface area parameters (because so far NET files with our naming convention only use these parameters for computing reaction parameters and initial species counts and we do not want to change reaction rates or concentrations)
parameters_names_backup_map_keys = parameters_names_backup_map.keys();
for i = 1:length(network_info.parameters)
    parameter_name = network_info_parameters_keys{i};
    parameter = network_info.parameters(parameter_name);
    for j = 1:length(parameters_names_backup_map_keys)
        parameters_names_backup_map_key = parameters_names_backup_map_keys{j};
        parameters_names_backup_map_value = parameters_names_backup_map(parameters_names_backup_map_key);
        parameter = parameter.replaceVariable(parameters_names_backup_map_key, parameters_names_backup_map_value);
    end
    updateParametersMap(parameter_name, parameter);
end
% Create eff_width parameter
updateParametersMap('eff_width', output_effective_width);

% updateParametersMap();
inferUnits();


for i = 1:length(network_info.reactions)
    reaction = network_info.reactions(i);
    for j = 1:length(parameters_names_backup_map_keys)
        parameters_names_backup_map_key = parameters_names_backup_map_keys{j};
        parameters_names_backup_map_value = parameters_names_backup_map(parameters_names_backup_map_key);
        reaction.rate_constant = reaction.rate_constant.replaceVariable(parameters_names_backup_map_key, parameters_names_backup_map_value);
    end
    network_info.reactions(i) = reaction;
end

for i = 1:length(network_info.species)
    single_species = network_info.species(i);
    for j = 1:length(parameters_names_backup_map_keys)
        parameters_names_backup_map_key = parameters_names_backup_map_keys{j};
        parameters_names_backup_map_value = parameters_names_backup_map(parameters_names_backup_map_key);
        single_species.diffusion_coefficient = single_species.diffusion_coefficient.replaceVariable(parameters_names_backup_map_key, parameters_names_backup_map_value);
        single_species.count = single_species.count.replaceVariable(parameters_names_backup_map_key, parameters_names_backup_map_value);
    end
    network_info.species(i) = single_species;
end

network_info_compartments_keys = network_info.compartments.keys;
for i = 1:length(network_info_compartments_keys)
    compartment_name = network_info_compartments_keys{i};
    compartment = network_info.compartments(compartment_name);
    for j = 1:length(parameters_names_backup_map_keys)
        parameters_names_backup_map_key = parameters_names_backup_map_keys{j};
        parameters_names_backup_map_value = parameters_names_backup_map(parameters_names_backup_map_key);
        compartment.size_expression = compartment.size_expression.replaceVariable(parameters_names_backup_map_key, parameters_names_backup_map_value);
    end
    network_info.compartments(compartment_name) = compartment;
end


% Add new volume parameters
for i = 1:length(all_compartment_data)
    object = all_compartment_data(i);
    name = object.name;
    object_volume = object.object_volume;
    units = output_volume_unit;
    % network_info_parameters([compartment_volume_prefix,name]) = DimensionedExpression(object_volume, units);
    updateParametersMap([compartment_volume_prefix,name], DimensionedExpression(object_volume, units));
end

% updateParametersMap();
inferUnits();


% Add new surface area parameters
for k = 1:length(all_membrane_data)
    object = all_membrane_data(k);
    name = object.name;
    area = object.area;
    area_parameter_name = [compartment_area_prefix,name];
    % network_info_parameters(area_parameter_name) = DimensionedExpression(area, output_area_unit);
    updateParametersMap(area_parameter_name, DimensionedExpression(area, output_area_unit));
    % network_info_parameters([compartment_volume_prefix,name]) = DimensionedExpression([area_parameter_name,'*eff_width'], output_volume_unit);
    updateParametersMap([compartment_volume_prefix,name], DimensionedExpression([area_parameter_name,'*eff_width'], output_volume_unit));
end

% updateParametersMap();
inferUnits();


%{
network_info_parameters_keys = network_info.parameters.keys';
% network_info_parameters_values = cell2mat(network_info.parameters.values');
network_info_parameters_values = cellfun(@(x)x, network_info.parameters.values);
network_info_parameters_values_index = [network_info_parameters_values.index]';
network_info_parameters_values_name = {network_info_parameters_values.name}';
network_info_parameters_values_value_expression = {network_info_parameters_values.value_expression}';
network_info_parameters_values_value_expression_is_numeric = cellfun(@(x)isStringNumeric(x), network_info_parameters_values_value_expression);
network_info_parameters_values_comment = {network_info_parameters_values.comment}';
network_info_parameters_values_units = {network_info_parameters_values.units}';

parameters_names_to_expressions_map = containers.Map(network_info_parameters_values_name, network_info_parameters_values_value_expression);
parameter_objects = struct('index', network_info_parameters_values_index, 'name', network_info_parameters_values_name, 'value', network_info_parameters_values_value_expression, 'is_numeric', num2cell(network_info_parameters_values_value_expression_is_numeric), 'comment', network_info_parameters_values_comment, 'units', network_info_parameters_values_units);
%}




% Set up VCML structure


warning('CellOrganizer:readNetworkIntoGeometry', 'Ignoring CSGdata.primitiveOnly.');
warning('CellOrganizer:readNetworkIntoGeometry', 'Ordinals not used here, connectivity inferred using adjacency in image.');


warning('CellOrganizer:readNetworkIntoGeometry', 'consensus_compartment: Check if this follows rules in Sekar and Faeder 2012');
function consensus_compartment = getConsensusCompartment(given_compartments)
    given_compartments = unique(given_compartments);
    if length(given_compartments) == 1
        consensus_compartment = given_compartments{1};
    else
        given_compartments_is_membrane = is_membrane_function(given_compartments);
        if length(given_compartments) == 2 && sum(given_compartments_is_membrane) == 0
            % Two compartments meet at the membrane
            consensus_compartment = get_compartment_pair_string(given_compartments{:}, name_map);
        elseif sum(given_compartments_is_membrane) == 1
            % Assume the VCML structure can be set to the membrane
            % TODO: Check that this is correct
            consensus_compartment = given_compartments{given_compartments_is_membrane};
        else
            consensus_compartment = given_compartments(given_compartments_is_membrane);
            consensus_compartment = consensus_compartment{1};
        end
    end
end


% Choose a compartment for each species
species_chosen_compartments = cell(length(network_info.species), 1);
for i = 1:length(network_info.species)
    single_species = network_info.species(i);
    % Get VCML structure name
    species_chosen_compartments{i} = getConsensusCompartment(single_species.compartments);
end
[network_info.species.chosen_compartment] = species_chosen_compartments{:};

% Compute initial concentration for each species
warning('CellOrganizer:readNetworkIntoGeometry', 'Assuming NET initial values are extensive as recommended for BNG');
species_initial_concentrations = repmat(DimensionedExpression(), 1, 0);
for i = 1:length(network_info.species)
    single_species = network_info.species(i);
    count = single_species.count;
    compartment_size_expression = network_info.compartments(single_species.chosen_compartment).size_expression;
    
    %{
    if isStringNumeric(compartment_size_expression)
        compartment_size_expression = str2double(compartment_size_expression);
    end
    %}
    
    % Don't check if this is a membrane or a volume
    % "Molecules in surfaces are assumed to be restricted to a small volume enveloping the surface (i.e., the surface volume is equal to the surface area multiplied by a surface thickness) provided by the modeler."
    
    %{
    if isnumeric(count)
        count = num2str(count);
    end
    if isnumeric(compartment_size_expression)
        compartment_size_expression = double2str(compartment_size_expression);
    end
    
    count = DimensionedExpression(count, 'molecule');
    compartment_size_expression = DimensionedExpression(compartment_size_expression, net_volume_unit);
    %}
    
    % Count to concentration
    concentration_conversion_factor = 1 / (parameterReference('N_Avo') * compartment_size_expression);
    concentration = count * concentration_conversion_factor;
    
    species_initial_concentrations(i) = concentration;
end
species_initial_concentrations = num2cell(species_initial_concentrations);
[network_info.species.concentration] = species_initial_concentrations{:};


% Add reaction names
reactions_digits = 0;
reactions_name_format = '';
function result = getReactionName(given_reaction)
    result = sprintf(reactions_name_format, given_reaction.index);
end
function result = getExtendedReactionName(given_reaction)
    result = strjoin({'reactants', strjoin(cellfun(@(x)num2str(x), num2cell(given_reaction.reactant_indices), 'UniformOutput', false), '_'), 'products', strjoin(cellfun(@(x)num2str(x), num2cell(given_reaction.product_indices), 'UniformOutput', false), '_')}, '_');
end
function updateReactionsNames()
    maximum_number_reactions = length(network_info.reactions);
    if add_translocation_intermediates
        % Include maximum number of intermediate reactions
        maximum_number_reactions = maximum_number_reactions * length(network_info.compartments);
    end
    reactions_digits = ceil(log10(maximum_number_reactions+1));
    reactions_name_format = ['r%0', num2str(reactions_digits), 'i'];
    reactions_names = arrayfun(@getReactionName, network_info.reactions, 'UniformOutput', false);
    reactions_extended_names = arrayfun(@getExtendedReactionName, network_info.reactions, 'UniformOutput', false);
    [network_info.reactions.name] = reactions_names{:};
    [network_info.reactions.extended_name] = reactions_extended_names{:};
end
updateReactionsNames();


% Decide on the compartment for each reaction
warning('CellOrganizer:readNetworkIntoGeometry', 'network_info.reactions.compartment values not known to be correct for BNG and/or VCell');
function result = getReactionConsensusCompartment(given_reaction)
    given_reaction_species_indices = [given_reaction.reactant_indices, given_reaction.product_indices];
    given_reaction_species = arrayfun(@(x)network_info.species(x), given_reaction_species_indices);
    given_reaction_species_compartments = cell(length(given_reaction_species), 1);
    for j = 1:length(given_reaction_species)
        given_reaction_species_compartments{j} = given_reaction_species(j).compartments;
    end
    given_reaction_species_compartments = cat(2, given_reaction_species_compartments{:});
    given_reaction_species_compartments = unique(given_reaction_species_compartments);
    result = getConsensusCompartment(given_reaction_species_compartments);
end
reactions_compartments = cell(length(network_info.reactions), 1);
for i = 1:length(network_info.reactions)
    reaction = network_info.reactions(i);
    reactions_compartments{i} = getReactionConsensusCompartment(reaction);
end
[network_info.reactions.compartment] = reactions_compartments{:};


% Decide on the reaction rate and units for each reaction
warning('CellOrganizer:readNetworkIntoGeometry', 'network_info.reactions.rate_constant_units values not known to be correct for BNG and/or VCell');
warning('CellOrganizer:readNetworkIntoGeometry', 'network_info.reactions.rate_constant values are not processed as microscopic');
reactions_rates = cell(length(network_info.reactions), 1);
reactions_rate_constants = cell(length(network_info.reactions), 1);
reactions_rates_units = cell(length(network_info.reactions), 1);
for i = 1:length(network_info.reactions)
    reaction = network_info.reactions(i);
    reaction_compartment = reaction.compartment;
    reaction_n_reactants = length(reaction.reactant_indices);
    reaction_rate_constant = reaction.rate_constant;
    
    % https://github.com/virtualcell/vcell/blob/master/vcell-core/src/main/java/cbit/vcell/model/GeneralKinetics.java#L109
    
    % https://github.com/virtualcell/vcell/blob/master/vcell-core/src/main/java/cbit/vcell/model/GeneralKinetics.java#L100
    % Sekar and Faeder, "Rule-Based Modeling of Signal Transduction: A Primer," in Computational Modeling of Signaling Networks
    % Correct volume from BNGL assuming the file follows guidelines
    % "If bimolecular and higher order reactions occur in different volumes, then the reaction rate constant needs to be scaled differently for the reaction occurring in each volume."
    % "If cBNG is invoked, i.e., the compartments and volumes are specified in the compartments block, then the volume scaling is part of the cBNG processing and the modeler needs to only provide the Avogadro number factor."
    % "In accordance with the explanation given in Subheading 5.1, unimolecular reactions are never scaled by volume, either in surfaces or volumes. Bimolecular reactions in 'volumes' and 'surfaces' are automatically scaled by the respective compartmental volume."
    % "It is not recommended to model reaction orders higher than bimolecular in the compartmental framework."
    % cBNGL-derived NET files appear to include these corrections
    %{
    if reaction_n_reactants > 1
        reaction_rate_constant = reaction_rate_constant / network_info.parameters([compartment_volume_prefix,reaction_compartment])^(reaction_n_reactants-1);
    elseif reaction_n_reactants < 1
        error('reaction_n_reactants = %0i', reaction_n_reactants);
    end
    reaction_rate_constant = reaction_rate_constant / DimensionedExpression(1, output_time_unit);
    %}
    
    reaction_rate = reaction_rate_constant;
    if reaction_n_reactants > 0
        for j = reaction.reactant_indices
            reaction_rate = reaction_rate * DimensionedExpression(network_info.species(j).name, output_concentration_unit);
        end
    end
    
    reactions_rates{i} = reaction_rate;
    reactions_rate_constants{i} = reaction_rate_constant;
end
[network_info.reactions.rate] = reactions_rates{:};
[network_info.reactions.rate_constant] = reactions_rate_constants{:};


% Split BNGL translocations because VCell does not accept reactions between non-adjacent compartments
reactions_to_add = network_info.reactions(1:0);
reactions_to_remove = [];
warning('CellOrganizer:readNetworkIntoGeometry', 'Splitting BNGL translocations with assumptions about diffusion coefficients and reaction rate, see code');
warning('CellOrganizer:readNetworkIntoGeometry', 'Permitting BNGL translocations that are invalid according to Sekar and Faeder 2012, "Rule-Based Modeling of Signal Transduction: A Primer" but used in the published example "journal.pcbi.1004611.s003.bngl" with assumptions about diffusion coefficients and reaction rate, see code');
for i = 1:length(network_info.reactions)
    reaction = network_info.reactions(i);
    
    reaction_reactant_species_indices = reaction.reactant_indices;
    reaction_product_species_indices = reaction.product_indices;
    reaction_species_indices = [reaction_reactant_species_indices, reaction_product_species_indices];
    
    reaction_reactant_species_compartments = cell(length(reaction_reactant_species_indices), 1);
    for j = 1:length(reaction_reactant_species_compartments)
        reaction_reactant_species_compartments{j} = {network_info.species(reaction_reactant_species_indices(j)).chosen_compartment};
    end
    reaction_reactant_species_compartments = cat(2, reaction_reactant_species_compartments{:});
    reaction_reactant_species_compartments = unique(reaction_reactant_species_compartments);
    
    reaction_product_species_compartments = cell(length(reaction_product_species_indices), 1);
    for j = 1:length(reaction_product_species_compartments)
        reaction_product_species_compartments{j} = {network_info.species(reaction_product_species_indices(j)).chosen_compartment};
    end
    reaction_product_species_compartments = cat(2, reaction_product_species_compartments{:});
    reaction_product_species_compartments = unique(reaction_product_species_compartments);
    
    reaction_species_compartments = [reaction_reactant_species_compartments, reaction_product_species_compartments];
    reaction_species_compartments = unique(reaction_species_compartments);
    % reaction_species_compartments
    if length(reaction_species_compartments) < 2
        % Do nothing
    elseif length(reaction_species_compartments) == 2
        % if is_compartment_pair_adjacent(reaction_species_compartments{1}, reaction_species_compartments{2}, 1)
        if is_network_info_compartment_pair_adjacent(reaction_species_compartments{1}, reaction_species_compartments{2}, 1)
            % Do nothing
        else
            if add_translocation_intermediates
                % Find path of compartments from reaction_species_compartments{1} to reaction_species_compartments{2}
                error('This section not updated to use UnitsManager and DimensionedExpression');
                compartment_chain1_connection_direction = '';
                compartment_chain1_insides = {reaction_species_compartments(1)};
                compartment_chain1_outside = {reaction_species_compartments{1}};
                compartment_chain2_insides = {reaction_species_compartments(2)};
                compartment_chain2_outside = {reaction_species_compartments{2}};
                while isempty(compartment_chain1_connection_direction)
                    if isempty(compartment_chain1_insides)
                        previous_compartment_chain1_insides = reaction_species_compartments(1);
                        previous_compartment_chain1_outside = reaction_species_compartments{1};
                        previous_compartment_chain2_insides = reaction_species_compartments(2);
                        previous_compartment_chain2_outside = reaction_species_compartments{2};
                    else
                        previous_compartment_chain1_insides = compartment_chain1_insides{end};
                        previous_compartment_chain1_outside = compartment_chain1_outside{end};
                        previous_compartment_chain2_insides = compartment_chain2_insides{end};
                        previous_compartment_chain2_outside = compartment_chain2_outside{end};
                    end
                    
                    if ~isempty(previous_compartment_chain1_insides)
                        current_compartment_chain1_insides = {};
                        for k = 1:length(previous_compartment_chain1_insides)
                            current_compartment_chain1_insides = [current_compartment_chain1_insides, getfield(network_info.compartments(previous_compartment_chain1_insides{k}), 'insides')];
                        end
                        compartment_chain1_insides{end+1} = current_compartment_chain1_insides;
                    end
                    
                    if ~isempty(previous_compartment_chain1_outside)
                        current_compartment_chain1_outside = network_info.compartments(previous_compartment_chain1_outside).outside;
                        if ~isempty(current_compartment_chain1_outside)
                            compartment_chain1_outside{end+1} = current_compartment_chain1_outside;
                        end
                    end
                    
                    if ~isempty(previous_compartment_chain2_insides)
                        current_compartment_chain2_insides = {};
                        for k = 1:length(previous_compartment_chain2_insides)
                            current_compartment_chain2_insides = [current_compartment_chain2_insides, getfield(network_info.compartments(previous_compartment_chain2_insides{k}), 'insides')];
                        end
                        compartment_chain2_insides{end+1} = current_compartment_chain2_insides;
                    end
                    
                    if ~isempty(previous_compartment_chain2_outside)
                        current_compartment_chain2_outside = network_info.compartments(previous_compartment_chain2_outside).outside;
                        if ~isempty(current_compartment_chain2_outside)
                            compartment_chain2_outside{end+1} = current_compartment_chain2_outside;
                        end
                    end
                    
                    insides_outside_matched = inCellArray(compartment_chain1_insides{end}, compartment_chain2_outside{end});
                    outside_insides_matched = inCellArray(compartment_chain2_insides{end}, compartment_chain1_outside{end});
                    if insides_outside_matched
                        compartment_chain1_connection_direction = 'insides';
                        break;
                    elseif outside_insides_matched
                        compartment_chain1_connection_direction = 'outside';
                        break;
                    elseif isempty(current_compartment_chain1_insides) && isempty(current_compartment_chain1_outside) && isempty(current_compartment_chain2_insides) && isempty(current_compartment_chain2_outside)
                        error('Path between non-adjacent compartments of reaction participants cannot be found');
                    end
                end
                
                if strcmp(compartment_chain1_connection_direction, 'insides')
                    compartment_chain = compartment_chain2_outside;
                    compartment_chain_end = reaction_species_compartments{1};
                else
                    compartment_chain = compartment_chain1_outside;
                    compartment_chain_end = reaction_species_compartments{2};
                end
                while ~isempty(compartment_chain{end}) && ~strcmp(compartment_chain{end}, compartment_chain_end)
                    compartment_chain_outside_compartment = network_info.compartments(compartment_chain{end}).outside;
                    if isempty(compartment_chain_outside_compartment)
                        break;
                    end
                    compartment_chain{end+1} = compartment_chain_outside_compartment;
                end
                
                
                % Build chain of reactions to approximate transport
                intermediate_species1_compartment_name = getConsensusCompartment(reaction_reactant_species_compartments);
                if strcmp(intermediate_species1_compartment_name, compartment_chain{1})
                    % Do nothing
                elseif strcmp(intermediate_species1_compartment_name, compartment_chain{end})
                    compartment_chain = compartment_chain(end:-1:1);
                else
                    error('intermediate_species1_compartment_name is not at either end of compartment_chain');
                end
                
                intermediate_species1 = struct;
                intermediate_species1.index = length(network_info.species) + 1;
                % intermediate_species1.species_graph = ['@',num3str(intermediate_species1.index)];
                intermediate_species1.species_graph = '';
                intermediate_species1.concentration = 0;
                intermediate_species1.compartments = {intermediate_species1_compartment_name};
                intermediate_species1.diffusion_coefficient = mean(getfield([network_info.species(reaction_species_indices)], 'diffusion_coefficient'));
                intermediate_species1.comment = ['intermediate species ',num2str(1)',' for translocation reaction "', reaction.name, '"'];
                intermediate_species1.name = species_index_to_name_function(intermediate_species1.index);
                intermediate_species1.extended_name = intermediate_species1.name;
                intermediate_species1.chosen_compartment = intermediate_species1_compartment_name;
                network_info.species(end+1) = intermediate_species1;
                
                intermediate_reaction1 = reaction;
                intermediate_reaction1.index = length(network_info.reactions) + length(reactions_to_add) + 1;
                intermediate_reaction1.product_indices = intermediate_species1.index;
                intermediate_reaction1.name = getReactionName(intermediate_reaction1);
                intermediate_reaction1.extended_name = getExtendedReactionName(intermediate_reaction1);
                intermediate_reaction1.compartment = getReactionConsensusCompartment(intermediate_reaction1);
                intermediate_reaction1.comment = ['intermediate reaction ',num2str(1)',' for translocation reaction "', reaction.name, '"'];
                reactions_to_add(end+1) = intermediate_reaction1;
                
                previous_intermediate_species2 = intermediate_species1;
                previous_intermediate_reaction2 = intermediate_reaction1;
                for k = 2:length(compartment_chain)
                    intermediate_species2_compartment_name = compartment_chain{k};
                    intermediate_species2 = struct;
                    intermediate_species2.index = length(network_info.species) + 1;
                    % intermediate_species2.species_graph = ['@',num3str(intermediate_species2.index)];
                    intermediate_species2.species_graph = '';
                    intermediate_species2.concentration = 0;
                    intermediate_species2.compartments = {intermediate_species2_compartment_name};
                    intermediate_species2.diffusion_coefficient = previous_intermediate_species2.diffusion_coefficient;
                    intermediate_species2.comment = ['intermediate species ',num2str(k)',' for translocation reaction "', reaction.name, '"'];
                    intermediate_species2.name = species_index_to_name_function(intermediate_species2.index);
                    intermediate_species2.extended_name = intermediate_species2.name;
                    intermediate_species2.chosen_compartment = intermediate_species2_compartment_name;
                    network_info.species(end+1) = intermediate_species2;
                    
                    intermediate_reaction2 = reaction;
                    intermediate_reaction2.index = length(network_info.reactions) + length(reactions_to_add) + 1;
                    intermediate_reaction2.reactant_indices = previous_intermediate_species2.index;
                    intermediate_reaction2.product_indices = intermediate_species2.index;
                    intermediate_reaction2.name = getReactionName(intermediate_reaction2);
                    intermediate_reaction2.extended_name = getExtendedReactionName(intermediate_reaction2);
                    intermediate_reaction2.compartment = getReactionConsensusCompartment(intermediate_reaction2);
                    intermediate_reaction2.comment = ['intermediate reaction ',num2str(k)',' for translocation reaction "', reaction.name, '"'];
                    
                    previous_intermediate_reaction2_rate = previous_intermediate_reaction2.rate;
                    intermediate_reaction2_rate = previous_intermediate_reaction2_rate;
                    % Replace reactants in reaction rate of previous_intermediate_reaction2
                    if isnumeric(intermediate_reaction2_rate)
                        intermediate_reaction2_rate = double2str(intermediate_reaction2_rate);
                    end
                    if length(previous_intermediate_reaction2.reactant_indices) ~= 1
                        intermediate_reaction2_rate = ['(',intermediate_reaction2_rate,')^(1/',num2str(length(previous_intermediate_reaction2.reactant_indices)),')'];
                    end
                    previous_intermediate_reaction2_product = network_info.species(previous_intermediate_reaction2.product_indices(1));
                    for k = 1:length(previous_intermediate_reaction2.reactant_indices)
                        previous_intermediate_reaction2_reactant = network_info.species(previous_intermediate_reaction2.reactant_indices(k));
                        intermediate_reaction2_rate = strrep(intermediate_reaction2_rate, previous_intermediate_reaction2_reactant.name, previous_intermediate_reaction2_product.name);
                    end
                    
                    % previous_intermediate_reaction2_rate
                    % intermediate_reaction2_rate
                    % error('Unfinished (add power of new reactant concentrations)')
                    
                    intermediate_reaction2.rate = intermediate_reaction2_rate;
                    % network_info.reactions(end+1) = intermediate_reaction2;
                    reactions_to_add(end+1) = intermediate_reaction2;
                    
                    % intermediate_reaction2
                    
                    previous_intermediate_species2 = intermediate_species2;
                    previous_intermediate_reaction2 = intermediate_reaction2;
                end
            end
            
            reactions_to_remove(end+1) = i;
        end
    else
        error('Reactions with reactants and products in more than two compartments not supported');
    end
end
reactions_to_keep = true(size(network_info.reactions));
reactions_to_keep(reactions_to_remove) = false;
network_info.reactions = network_info.reactions(reactions_to_keep);
network_info.reactions(end+1:end+length(reactions_to_add)) = reactions_to_add;
% updateSpeciesNames();
% updateReactionsNames();


% reactions_rate_units = {network_info.reactions.rate_units}';

warning('CellOrganizer:readNetworkIntoGeometry', 'Rate constant units converted improperly');
%{
species_rates = repmat(DimensionedExpression(0), length(network_info.species), 1);
reaction_mass_action_rates = repmat(DimensionedExpression(0), length(network_info.species), 1);
for i = 1:length(network_info.reactions)
    reaction = network_info.reactions(i);
    
    % Reactant nodes and stoichiometries
    reaction_reactant_indices_unique = unique(reaction.reactant_indices);
    reactant_stoichiometries = arrayfun(@(x)sum(reaction.reactant_indices == x), reaction_reactant_indices_unique);
    reactant_first_indices = arrayfun(@(x)find(reaction.reactant_indices == x, 1, 'first'), reaction_reactant_indices_unique);
    
    % Product nodes and stoichiometries
    reaction_product_indices_unique = unique(reaction.product_indices);
    product_stoichiometries = arrayfun(@(x)sum(reaction.product_indices == x), reaction_product_indices_unique);
    product_first_indices = arrayfun(@(x)find(reaction.product_indices == x, 1, 'first'), reaction_product_indices_unique);
    
    % Construct species rate expressions
    for j = 1:length(reaction_reactant_indices_unique)
        reactant_index = reaction.reactant_indices(reactant_first_indices(j));
        species_rates(reactant_index) = species_rates(reactant_index) - reactant_stoichiometries(j) * reaction.rate;
    end
    for j = 1:length(reaction_product_indices_unique)
        product_index = reaction.product_indices(product_first_indices(j));
        species_rates(product_index) = species_rates(product_index) + product_stoichiometries(j) * reaction.rate;
    end
    
    % Construct reaction rate expression for mass action kinetics
    reaction_mass_action_rate = reaction_mass_action_rates(i);
    reaction_mass_action_rate_forward = DimensionedExpression('Kf');
    for j = 1:length(reaction_reactant_indices_unique)
        reactant_index = reaction.reactant_indices(reactant_first_indices(j));
        reaction_mass_action_rate_forward = reaction_mass_action_rate_forward * network_info.species(reactant_index).name;
    end
    reaction_mass_action_rate_reverse = DimensionedExpression('Kr');
    for j = 1:length(reaction_product_indices_unique)
        product_index = reaction.product_indices(product_first_indices(j));
        reaction_mass_action_rate_reverse = reaction_mass_action_rate_reverse * network_info.species(product_index).name;
    end
    reaction_mass_action_rate = reaction_mass_action_rate_forward - reaction_mass_action_rate_reverse;
    reaction_mass_action_rates(i) = reaction_mass_action_rate;
end
warning('CellOrganizer:readNetworkIntoGeometry', 'ParameterNode with role ''user defined'' attribute Unit not computed, assuming ??? for 1 parameter');
species_rates = num2cell(species_rates);
[network_info.species.rate] = species_rates{:};
reaction_mass_action_rates = num2cell(reaction_mass_action_rates);
[network_info.reactions.mass_action_rate] = reaction_mass_action_rates{:};

% species_rates
% error('Unfinished')
%}

warning('CellOrganizer:readNetworkIntoGeometry', 'Unit conversions unfinished');


% Add MembraneVoltage constants here so they are in MathDescription Constants but not in ModelParameters Parameters
for j = 1:length(all_membrane_data)
    object = all_membrane_data(j);
    updateParametersMap(object.MembraneVoltageName, DimensionedExpression(0, 'V'), false);
end

% updateParametersMap();
inferUnits();


network_info.parameters_names_topological_order = parameters_names_topological_order;
network_info.parameters_should_write_to_parameters = parameters_should_write_to_parameters;
network_info.set_should_write_to_parameters = @set_should_write_to_parameters;
network_info.get_should_write_to_parameters = @get_should_write_to_parameters;



warning('CellOrganizer:readNetworkIntoGeometry', 'TODO: Double check unit adjustments (um_to_output_length_unit, etc.)');


geometry_info = struct();
geometry_info.models = models;
geometry_info.imgs = imgs;
geometry_info.meshes = {};
geometry_info.meshData = meshData;
geometry_info.dim_chars = dim_chars;
geometry_info.sign_chars = sign_chars;
geometry_info.connectivity = connectivity;
geometry_info.connectivity_se = connectivity_se;
geometry_info.avogadro_constant = avogadro_constant;
geometry_info.avogadro_constant_value_expression = avogadro_constant_value_expression;
geometry_info.avogadro_constant_units = avogadro_constant_units;
geometry_info.get_compartment_pair_string = @get_compartment_pair_string;
geometry_info.separate_compartment_images = @separate_compartment_images;
geometry_info.combine_compartment_images = @combine_compartment_images;
geometry_info.get_compartment_image_boundary_value = @get_compartment_image_boundary_value;
geometry_info.crop_compartment_image = @crop_compartment_image;
geometry_info.combine_resize_compartment_images = @combine_resize_compartment_images;
geometry_info.model_names = model_names;
geometry_info.model_cytonuclearflags = model_cytonuclearflags;
geometry_info.model_names_translated = model_names_translated;
geometry_info.name_map = name_map;
geometry_info.framework_compartment_names = framework_compartment_names;
geometry_info.is_framework_compartment_name_function = is_framework_compartment_name_function;
geometry_info.named_imgs = named_imgs;
geometry_info.named_imgs_volumes = named_imgs_volumes;
geometry_info.named_imgs_exclusive_volumes = named_imgs_exclusive_volumes;
geometry_info.named_imgs_cytonuclearflags = named_imgs_cytonuclearflags;
geometry_info.names = names;
% geometry_info.named_imgs_volumes_keys = named_imgs_volumes_keys;
% geometry_info.named_imgs_volumes_keys_sorted = named_imgs_volumes_keys_sorted;
geometry_info.all_compartments_image_before_downsampling = all_compartments_image_before_downsampling;
geometry_info.all_compartments_image = all_compartments_image;
geometry_info.named_imgs_before_downsampling = named_imgs_before_downsampling;
% geometry_info.resolution_before_downsampling = resolution_before_downsampling;
geometry_info.resolution_single_before_downsampling = resolution_single_before_downsampling;
geometry_info.resolution = resolution;
geometry_info.resolution_single = resolution_single;
% geometry_info.adjacent_pairs_before_downsampling = adjacent_pairs_before_downsampling;
geometry_info.adjacent_pairs = adjacent_pairs;
geometry_info.names_sorted = names_sorted;
geometry_info.n_compartments = n_compartments;
geometry_info.all_compartment_volumes = all_compartment_volumes;
geometry_info.all_compartment_exclusive_volumes = all_compartment_exclusive_volumes;
geometry_info.all_compartment_indices = all_compartment_indices;
geometry_info.n_network_info_compartments = n_network_info_compartments;
geometry_info.network_info_compartments_keys = network_info_compartments_keys;
geometry_info.network_info_compartments_keys_indices_map = network_info_compartments_keys_indices_map;
geometry_info.network_info_adjacency_max_degree = network_info_adjacency_max_degree;
geometry_info.network_info_adjacency_matrix = network_info_adjacency_matrix;
geometry_info.network_info_adjacency_matrices = network_info_adjacency_matrices;
geometry_info.is_network_info_compartment_pair_adjacent = @is_network_info_compartment_pair_adjacent;
geometry_info.adjacent_pairs_max_degree = adjacent_pairs_max_degree;
geometry_info.adjacent_pairs_matrix = adjacent_pairs_matrix;
geometry_info.adjacent_pairs_matrices = adjacent_pairs_matrices;
geometry_info.is_compartment_pair_adjacent = @is_compartment_pair_adjacent;
geometry_info.all_compartments_objects = all_compartments_objects;
geometry_info.all_compartments_objects_image = all_compartments_objects_image;
geometry_info.all_compartments_objects_names = all_compartments_objects_names;
geometry_info.all_compartments_objects_names_to_compartment_names = all_compartments_objects_names_to_compartment_names;
geometry_info.all_compartments_objects_indices = all_compartments_objects_indices;
geometry_info.n_compartments_objects = n_compartments_objects;
geometry_info.object_adjacent_pairs = object_adjacent_pairs;
geometry_info.getAdjacentValues = @getAdjacentValues;
geometry_info.getAdjacentCompartments = @getAdjacentCompartments;

geometry_info.all_compartment_data = all_compartment_data;
geometry_info.all_object_data = all_object_data;
geometry_info.all_object_membrane_data = all_object_membrane_data;
geometry_info.all_membrane_data = all_membrane_data;
geometry_info.all_membrane_names = all_membrane_names;
geometry_info.is_membrane_function = @is_membrane_function;
% geometry_info.parameter_objects = parameter_objects;
geometry_info.species_index_to_name_function = @species_index_to_name_function;
% geometry_info.updateSpeciesNames = @updateSpeciesNames;
% geometry_info.expression_evaluation_function = @expression_evaluation_function;

geometry_info.translateWithDefaultIdentity = @translateWithDefaultIdentity;
geometry_info.parameterReferenceBase = @parameterReferenceBase;
geometry_info.evaluateExpression = @evaluateExpression;


network_with_geometry_info = struct();
network_with_geometry_info.network_info = network_info;
network_with_geometry_info.geometry_info = geometry_info;

end




function adjacent_pairs = getAdjacentValues(image, connectivity)
% Find adjacent values in image in directions

if nargin < 2
    connectivity = 26;
end
if connectivity == 6
    directions = [1, 0, 0; 0, 1, 0; 0, 0, 1];
elseif connectivity == 18
    directions = [1, 0, 0; 0, 1, 0; 0, 0, 1; 1, 1, 0; 1, -1, 0; 1, 0, 1; 1, 0, -1; 0, 1, 1; 0, 1, -1];
elseif connectivity == 26
    directions = [1, 0, 0; 0, 1, 0; 0, 0, 1; 1, 1, 0; 1, -1, 0; 1, 0, 1; 1, 0, -1; 0, 1, 1; 0, 1, -1; 1, 1, 1; 1, 1, -1; 1, -1, 1; -1, 1, 1];
else
    error('''directions'' must be 6, 18, or 26');
end
adjacent_pairs = {};
for direction = directions'
    temp = image(1 - min(direction(1), 0):end - max(direction(1), 0), 1 - min(direction(2), 0):end - max(direction(2), 0), 1 - min(direction(3), 0):end - max(direction(3), 0));
    temp2 = image(1 + max(direction(1), 0):end + min(direction(1), 0), 1 + max(direction(2), 0):end + min(direction(2), 0), 1 + max(direction(3), 0):end + min(direction(3), 0));
    temp = [temp(:), temp2(:)];
    temp = unique(temp, 'rows');
    clear temp2;
    adjacent_pairs{end+1} = temp;
end
adjacent_pairs = cell2mat(adjacent_pairs');
adjacent_pairs = sortrows(adjacent_pairs);
adjacent_pairs = sort(adjacent_pairs, 2);
adjacent_pairs = unique(adjacent_pairs, 'rows');
adjacent_pairs = adjacent_pairs(adjacent_pairs(:, 1) ~= adjacent_pairs(:, 2), :);

end




function adjacent_pairs = getAdjacentCompartments(network_info_compartments, all_compartment_indices)
% Find adjacent values for volumetric compartments using NET file info

adjacent_pairs = zeros(0, 2);
network_info_compartments_keys = network_info_compartments.keys();
all_compartment_indices_double_insides = containers.Map('KeyType', 'char', 'ValueType', 'any');
all_compartment_indices_double_outsides = containers.Map('KeyType', 'char', 'ValueType', 'any');
for i = 1:length(network_info_compartments)
    all_compartment_indices_double_insides(network_info_compartments_keys{i}) = {};
    all_compartment_indices_double_outsides(network_info_compartments_keys{i}) = {};
end
for i = 1:length(network_info_compartments)
    compartment = network_info_compartments(network_info_compartments_keys{i});
    compartment_name = compartment.name;
    for j = 1:length(compartment.insides)
        compartment_inside = compartment.insides{j};
        all_compartment_indices_double_insides(compartment_name) = [all_compartment_indices_double_insides(compartment_name), network_info_compartments(compartment_inside).insides];
    end
    all_compartment_indices_double_insides(compartment_name) = unique(all_compartment_indices_double_insides(compartment_name));
    if ~isempty(compartment.outside)
        all_compartment_indices_double_outsides(compartment_name) = network_info_compartments(compartment.outside).outside;
    end
end
for i = 1:length(network_info_compartments)
    compartment = network_info_compartments(network_info_compartments_keys{i});
    if compartment.spatial_dimensions == 3 && ~isempty(all_compartment_indices_double_outsides(compartment.name))
        adjacent_pairs(i, :) = [all_compartment_indices(compartment.name), all_compartment_indices(all_compartment_indices_double_outsides(compartment.name))];
    end
end
adjacent_pairs = sortrows(adjacent_pairs);
adjacent_pairs = sort(adjacent_pairs, 2);
adjacent_pairs = unique(adjacent_pairs, 'rows');
adjacent_pairs = adjacent_pairs(adjacent_pairs(:, 1) ~= adjacent_pairs(:, 2), :);

end




function adjacent_pairs = getAdjacentObjects(objects)
% Find adjacent objects for volumetric compartments using NET file, CSG, and mesh info

adjacent_pairs = zeros(0, 2);
objects_keys = objects.keys();
objects_indices_to_names = {};
for i = 1:length(objects)
    object = objects(objects_keys{i});
    object_index = object.index;
    object_outside = object.outside;
    object_insides = object.insides;
    objects_indices_to_names{object_index} = object.name;
    if ~isempty(object_outside)
        adjacent_pairs(end+1, :) = [object_index, objects(object_outside).index];
    end
    for j = 1:length(object_insides)
        object_inside = object_insides{j};
        adjacent_pairs(end+1, :) = [object_index, objects(object_inside).index];
    end
end
adjacent_pairs = sortrows(adjacent_pairs);
adjacent_pairs = sort(adjacent_pairs, 2);
adjacent_pairs = unique(adjacent_pairs, 'rows');
adjacent_pairs = adjacent_pairs(adjacent_pairs(:, 1) ~= adjacent_pairs(:, 2), :);

end




function value = inCellArray(cell_array, key)

value = any(strcmp(cell_array, key));

end




function value = intersectCellArrays(cell_array1, cell_array2)

value = {};
for i = 1:length(cell_array1)
    if inCellArray(cell_array1{i}, cell_array2)
        value{end+1} = cell_array1{i};
    end
end

end

