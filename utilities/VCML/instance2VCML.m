function [ result ] = instance2VCML( models, imgs, options, savepath )
%INSTANCE2VCML Writes geometry and compartmental reactions to VCell VCML file.
%
% Inputs
% ------
% models   = cell array of models produced in slml2img.m
% imgs     = cell array containing image arrays for all objects to be saved
% options  = options struct produced in slml2img.m
% savepath = output filename
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
% Copyright (C) 2012-2019 Murphy Lab
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

% 2018-10-02 Taraz Buck - Copied `instance2SBML3_mod.m` to `instance2VCML.m`.


debug = false;
% debug = true;

if ~debug
    warning('off', 'CellOrganizer:instance2VCML');
end

result = false;

if nargin < 3
    warning('CellOrganizer:instance2VCML', 'First three arguments are required');
    return;
end

if nargin < 4
    warning('CellOrganizer:instance2VCML', 'Argument savepath not given, defaulting to ''./model.vcml''');
    savepath = './model.vcml';
    [~, savepath_file, ~] = fileparts(savepath);
end

imgs_valid = true(length(imgs), 1);
for i = 1:length(imgs)
    imgs_valid(i) = imgs_valid(i) && numel(imgs{i}) > 0;
    imgs_valid(i) = imgs_valid(i) && ndims(imgs{i}) == 3;
end
if ~all(imgs_valid)
    warning('CellOrganizer:instance2VCML', 'Argument imgs must contain nonempty 3D arrays');
    return;
end

NETfile = options.output.NET.filename;
VCMLfile = options.output.VCML.input_filename;

VCMLTranslations = options.output.VCML.translations;
VCMLDownsample = options.output.VCML.downsampling;
VCMLAddTranslocationIntermediates = options.output.VCML.addTranslocationIntermediates;
% VCMLObjectsAlwaysPresent = options.output.VCML.objectsAlwaysPresent;

VCMLNumSimulations = 1;
% VCMLNumSimulations = options.output.VCML.numSimulations;
VCMLEndTime = options.output.VCML.endTime;
VCMLDefaultTimeStep = options.output.VCML.defaultTimeStep;
VCMLMinTimeStep = options.output.VCML.minTimeStep;
VCMLMaxTimeStep = options.output.VCML.maxTimeStep;
VCMLOutputTimeStep = options.output.VCML.outputTimeStep;
VCMLAbsoluteError = options.output.VCML.absoluteError;
VCMLRelativeError = options.output.VCML.relativeError;
VCMLDefaultDiffusionCoefficient = options.output.VCML.defaultDiffusionCoefficient;

net_length_unit = options.output.NET.units.length;
net_area_unit = [net_length_unit,'2'];
net_volume_unit = [net_length_unit,'3'];
net_time_unit = options.output.NET.units.time;
net_count_unit = 'molecules';
net_concentration_unit = options.output.NET.units.concentration;
warning('CellOrganizer:instance2VCML', 'Assuming NET file concentrations are extensive');

avogadro_constant = 6.022140857e23;

vcml_length_unit = 'um';
vcml_area_unit = [vcml_length_unit,'2'];
vcml_volume_unit = [vcml_length_unit,'3'];
vcml_concentration_unit = 'uM';
vcml_volume_substance_unit = 'uM.um3';
vcml_membrane_substance_unit = 'molecules';
vcml_time_unit = 's';
vcml_volume_reaction_rate_unit = strjoin({vcml_concentration_unit, [vcml_time_unit,'-1']},'.');
vcml_membrane_reaction_rate_unit = strjoin({vcml_membrane_substance_unit, [vcml_length_unit,'-2'], [vcml_time_unit,'-1']},'.');
vcml_diffusion_coefficient_unit = 'um2.s-1';

warning('CellOrganizer:instance2VCML', 'Assuming options.resolution is in um');
resolution = unit_convert('um', vcml_length_unit, options.resolution.cubic);

SI_effective_width = options.output.NET.effectiveWidth;
vcml_effective_width = unit_convert('m', vcml_length_unit, SI_effective_width);

% use_image_adjacency = options.output.NET.useImageAdjacency;
use_image_adjacency = false;



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

% kinetics_type = 'GeneralKinetics';
kinetics_type = 'MassAction';

avogadro_constant_value_expression = '6.022140857e23';
avogadro_constant_units = 'molecules.mol-1';

vcml_diffusion_coefficient = unit_convert('m2.s-1', vcml_diffusion_coefficient_unit, VCMLDefaultDiffusionCoefficient);

is_named_imgs_inclusive = true;

enforce_parent_buffer = true;
enforce_sibling_buffer = true;




% Process BioNetGen NET file

network_info = readNetwork(NETfile);


% Helper functions

function x = isStringNumeric(x)
    x = ~isnan(str2double(x));
end


function [given_pair_string_translated, given_pair_string_cell, given_pair_string] = get_compartment_pair_string(given_name1, given_name2, given_name_map)
    given_pair_string_cell = {translateWithDefaultIdentity(given_name_map, given_name1), translateWithDefaultIdentity(given_name_map, given_name2)};
    given_pair_string_cell = sort(given_pair_string_cell);
    given_pair_string = [given_pair_string_cell{1}, '_', given_pair_string_cell{2}];
    given_pair_string_translated = translateWithDefaultIdentity(given_name_map, given_pair_string);
end


function [given_named_imgs, given_named_imgs_volumes, given_names_sorted] = separate_compartment_images(given_all_compartments_image, given_resolution, given_names_sorted, given_network_info_compartments)
    % Recreate named_imgs by separating indexed regions in given_all_compartments_image and filling blank portions of parents with child regions
    
    % Separate indexed regions
    given_named_imgs = containers.Map('KeyType', 'char', 'ValueType', 'any');
    given_named_imgs_volumes = containers.Map('KeyType', 'char', 'ValueType', 'double');
    for gi = 1:length(given_names_sorted)
        given_name = given_names_sorted{gi};
        given_img = given_all_compartments_image == gi;
        given_named_imgs(given_name) = given_img;
        given_named_imgs_volumes(given_name) = sum(given_img(:)) * prod(given_resolution);
    end
    
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


function [given_all_compartments_image, given_named_imgs, given_names_sorted] = combine_compartment_images(given_named_imgs, given_names_sorted, given_resolution, given_network_info_compartments)
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
    
    [given_named_imgs, given_named_imgs_volumes, given_names_sorted] = separate_compartment_images(given_all_compartments_image, given_resolution, given_names_sorted, given_network_info_compartments);
    
    
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
        
        [given_named_imgs, given_named_imgs_volumes, given_names_sorted] = separate_compartment_images(given_all_compartments_image, given_resolution, given_names_sorted, given_network_info_compartments);
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


function [given_all_compartments_image, given_named_imgs, given_resolution, given_named_imgs_volumes, given_names_sorted] = combine_resize_compartment_images(given_named_imgs, given_resolution, given_names_sorted, given_network_info_compartments, given_scale)
    % Combine named_imgs into an indexed image, resample, and reproduce named_imgs
    
    [given_all_compartments_image, given_named_imgs, given_names_sorted] = combine_compartment_images(given_named_imgs, given_names_sorted, given_resolution, given_network_info_compartments);
    
    given_scale_pad_size = ceil(1 ./ given_scale);
    
    given_all_compartments_image = padarray(given_all_compartments_image, given_scale_pad_size, 'replicate');
    given_all_compartments_image_before_downsampling_size = size(given_all_compartments_image);
    given_all_compartments_image = image_resize_nd(given_all_compartments_image, given_scale, 'nearest');
    given_all_compartments_image_size = size(given_all_compartments_image);
    
    given_actual_scale = given_all_compartments_image_size ./ given_all_compartments_image_before_downsampling_size;
    given_resolution = given_resolution ./ given_actual_scale;
    
    % Assumes boundary should be one value and is of arbitrary size
    given_all_compartments_image = crop_compartment_image(given_all_compartments_image);
    
    [given_named_imgs, given_named_imgs_volumes, given_names_sorted] = separate_compartment_images(given_all_compartments_image, given_resolution, given_names_sorted, given_network_info_compartments);
    [given_all_compartments_image, given_named_imgs, given_names_sorted] = combine_compartment_images(given_named_imgs, given_names_sorted, given_resolution, given_network_info_compartments);
    [given_named_imgs, given_named_imgs_volumes, given_names_sorted] = separate_compartment_images(given_all_compartments_image, given_resolution, given_names_sorted, given_network_info_compartments);
end



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



% Create a map for translating compartment names

name_map = containers.Map('KeyType', 'char', 'ValueType', 'char');
name_map('EC') = 'EC';
name_map('cell') = 'cell';
name_map('nucleus') = 'nucleus';
% Defaults
for i = 1:length(models)
    name_map(model_names{i}) = model_names{i};
end
% Translations in options
for i = 1:size(VCMLTranslations, 1)
    name_map(VCMLTranslations{i, 1}) = VCMLTranslations{i, 2};
end


EC_translated = translateWithDefaultIdentity(name_map, 'EC');
cell_translated = translateWithDefaultIdentity(name_map, 'cell');
nucleus_translated = translateWithDefaultIdentity(name_map, 'nucleus');
EC_cell_translated = get_compartment_pair_string('EC', 'cell', name_map);
cell_nucleus_translated = get_compartment_pair_string('cell', 'nucleus', name_map);
framework_compartment_names = {EC_translated, cell_translated, nucleus_translated};
is_framework_compartment_name_function = @(x)any(strcmp(x, framework_compartment_names));

model_names_translated = cell(length(models), 1);
for i = 1:length(model_names)
    model_names_translated{i} = translateWithDefaultIdentity(name_map, model_names{i});
end


% Create a single image containing all compartments. Assumes no overlapping.

% Collect names and properties of compartments
warning('CellOrganizer:instance2VCML', 'Volumes might be handled differently than in BioNetGen and other software; compare https://github.com/RuleWorld/bionetgen/blob/master/bng2/Perl2/Compartment.pm');
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


% Add missing compartments to network_info.compartments
if any(strcmpi(options.synthesis, {'all', 'framework', 'cell'}))
    if ~network_info.compartments.isKey(EC_translated)
        network_info.compartments(EC_translated) = struct('name', EC_translated, 'spatial_dimensions', 3, 'size_expression', nan, 'outside', '', 'insides', {{EC_cell_translated}});
    end
    if ~network_info.compartments.isKey(EC_cell_translated)
        network_info.compartments(EC_cell_translated) = struct('name', EC_cell_translated, 'spatial_dimensions', 2, 'size_expression', nan, 'outside', EC_translated, 'insides', {{cell_translated}});
    end
    if ~network_info.compartments.isKey(cell_translated)
        temp = struct('name', cell_translated, 'spatial_dimensions', 3, 'size_expression', nan, 'outside', EC_cell_translated, 'insides', {{}});
        if any(strcmpi(options.synthesis, {'all', 'framework'}))
            temp.insides = {cell_nucleus_translated};
        end
        network_info.compartments(cell_translated) = temp;
    end
end
if ~network_info.compartments.isKey(nucleus_translated) && any(strcmpi(options.synthesis, {'all', 'framework', 'nucleus'}))
    temp = struct('name', nucleus_translated, 'spatial_dimensions', 3, 'size_expression', nan, 'outside', '', 'insides', {{}});
    if any(strcmpi(options.synthesis, {'all', 'framework'}))
        temp.outside = cell_nucleus_translated;
    end
    network_info.compartments(nucleus_translated) = temp;
end
if ~network_info.compartments.isKey(cell_nucleus_translated) && any(strcmpi(options.synthesis, {'all', 'framework'}))
    network_info.compartments(cell_nucleus_translated) = struct('name', cell_nucleus_translated, 'spatial_dimensions', 2, 'size_expression', nan, 'outside', cell_translated, 'insides', {{nucleus_translated}});
end


if any(strcmpi(options.synthesis, {'all', 'cell'}))
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
        network_info.compartments(name) = struct('name', name, 'spatial_dimensions', 3, 'size_expression', nan, 'outside', outside, 'insides', {{}});
        network_info.compartments(outside) = struct('name', outside, 'spatial_dimensions', 2, 'size_expression', nan, 'outside', outside_outside, 'insides', {{name}});
        compartment_outside_outside = network_info.compartments(outside_outside);
        compartment_outside_outside.insides{end+1} = outside;
        network_info.compartments(outside_outside) = compartment_outside_outside;
    end
end


[all_compartments_image, named_imgs, names_sorted] = combine_compartment_images(named_imgs, names_sorted, resolution, network_info.compartments);


for i = 1:length(names)
    name = names{i};
    img = named_imgs(name);
    named_imgs_volumes(name) = sum(img(:)) * prod(resolution);
end


all_compartment_names = names_sorted;
n_compartments = length(all_compartment_names);
% Exclusive volumes (do not include volumes of compartments inside)
all_compartment_volumes = zeros(n_compartments, 1);
for i = 2:n_compartments
    all_compartment_volumes(i) = named_imgs_volumes(all_compartment_names{i});
end

all_compartment_indices = containers.Map('KeyType', 'char', 'ValueType', 'double');
for i = 1:n_compartments
    all_compartment_indices(all_compartment_names{i}) = i;
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


% Create indexed image

% Combine first
[all_compartments_image, named_imgs, names_sorted] = combine_compartment_images(named_imgs, names_sorted, resolution, network_info.compartments);

if any(VCMLDownsample ~= 1)
    all_compartments_image_before_downsampling = all_compartments_image;
    named_imgs_before_downsampling = named_imgs;
    resolution_before_downsampling = resolution;
    
    [all_compartments_image, named_imgs, resolution, named_imgs_volumes, names_sorted] = combine_resize_compartment_images(named_imgs, resolution, names_sorted, network_info.compartments, ones(1, 3) .* VCMLDownsample);
end

if options.output.indexedimage
    if any(VCMLDownsample ~= 1)
        imwrite(reshape_contrast(single(all_compartments_image_before_downsampling), -1), [savepath, ' all_compartments_image_before_downsampling.png']);
    end
    imwrite(reshape_contrast(single(all_compartments_image), -1), [savepath, ' all_compartments_image.png']);
end


% Assign volumes and surface areas in network_info.compartments
for i = 1:length(named_imgs)
    name = names_sorted{i};
    compartment = network_info.compartments(name);
    compartment.size_expression = named_imgs_volumes(name);
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
fprintf('Compartment volumes:\n');
for i = 1:n_compartments
    fprintf('    %s: %.6f %s\n', all_compartment_names{i}, named_imgs_volumes(all_compartment_names{i}), vcml_volume_unit);
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

function result2 = is_network_info_compartment_pair_adjacent(given_compartment1, given_compartment2, given_degrees)
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
    result2 = network_info_adjacency_matrices(given_compartment1_index, given_compartment2_index, given_degrees);
end


% Create compartment adjacency matrices for paths of length 1 to adjacent_pairs_max_degree
adjacent_pairs_max_degree = 2;
adjacent_pairs_matrix = false(n_compartments, n_compartments);
adjacent_pairs_matrices = repmat(adjacent_pairs_matrix, 1, 1, adjacent_pairs_max_degree);
for i = 2:adjacent_pairs_max_degree
    adjacent_pairs_matrices(:, :, i) = single(adjacent_pairs_matrices(:, :, i-1)) * adjacent_pairs_matrix > 0;
end
function result2 = is_compartment_pair_adjacent(given_compartment1, given_compartment2, given_degrees)
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
    result2 = adjacent_pairs_matrices(given_compartment1_index, given_compartment2_index, given_degrees);
end


% Find all objects in each compartment

all_compartments_objects = containers.Map('KeyType', 'char', 'ValueType', 'any');
% New indexed image for individual compartmental objects
all_compartments_objects_image = [];
all_compartments_objects_names = {};
all_compartments_objects_names_to_compartment_names = containers.Map('KeyType', 'char', 'ValueType', 'char');
all_compartments_objects_indices = containers.Map('KeyType', 'char', 'ValueType', 'double');
should_only_keep_largest_framework_objects = true;
k = 1;
for i = 1:n_compartments
    name = all_compartment_names{i};
    name_img = named_imgs(name);
    name_img_objects = bwconncomp(name_img, connectivity);
    name_img_objects_stats = regionprops(name_img_objects);
    
    should_only_keep_largest_compartment_objects = should_only_keep_largest_framework_objects && is_framework_compartment_name_function(name) && length(name_img_objects_stats) > 1;
    if should_only_keep_largest_compartment_objects
        compartment_largest_object_index = nan;
        compartment_largest_object_volume = 0;
    end
    
    previous_all_compartments_objects_image = all_compartments_objects_image;
    for j = 1:length(name_img_objects_stats)
        compartment_object_name = [name, num2str(k)];
        
        mask = false(size(all_compartments_image));
        mask(name_img_objects.PixelIdxList{j}) = true;
        mask = mask & (all_compartments_image == all_compartment_indices(name));
        
        temp = name_img_objects_stats(j);
        temp.SamplePixel = name_img_objects.PixelIdxList{j}(1);
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
        
        if temp.FilledVolume == numel(mask)
            temp.FilledVolumeUnlessAll = temp.Volume;
            temp.FilledAreaUnlessAll = temp.Area;
        else
            temp.FilledVolumeUnlessAll = temp.FilledVolume;
            temp.FilledAreaUnlessAll = temp.FilledArea;
        end
        
        temp.Volume = temp.Volume * prod(resolution);
        temp.FilledVolume = temp.FilledVolume * prod(resolution);
        temp.FilledVolumeUnlessAll = temp.FilledVolumeUnlessAll * prod(resolution);
        % TODO: Correctly compute area of scaled geometry instead of approximating it by scaling voxel counts likes this
        temp.Area = temp.Area * power(prod(resolution), 2/3);
        temp.FilledArea = temp.FilledArea * power(prod(resolution), 2/3);
        temp.FilledAreaUnlessAll = temp.FilledAreaUnlessAll * power(prod(resolution), 2/3);
        
        temp.compartment_name = name;
        temp.compartment_object_name = compartment_object_name;
        temp.object_index = k;
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
            compartment_largest_object_index = temp.object_index;
            compartment_largest_object_volume = temp.Volume;
        end
        
        if ~should_only_keep_largest_compartment_objects || compartment_largest_object_index == temp.object_index
            all_compartments_objects(compartment_object_name) = temp;
            
            % Image
            if isempty(all_compartments_objects_image)
                all_compartments_objects_image = zeros(size(all_compartments_image), 'uint32');
            end
            
            all_compartments_objects_image(mask) = k;
            all_compartments_objects_names{k} = compartment_object_name;
            all_compartments_objects_names_to_compartment_names(compartment_object_name) = name;
            all_compartments_objects_indices(compartment_object_name) = k;
            
            k = k + 1;
        end
    end
end
n_compartments_objects = length(all_compartments_objects_names);


% Find adjacent compartment objects
if use_image_adjacency
    % Find adjacent compartments
    object_adjacent_pairs = getAdjacentValues(all_compartments_objects_image, connectivity);
else
    if ~should_only_keep_largest_framework_objects
        error('Not implemented');
    end
    object_adjacent_pairs = zeros(0, 2);
    all_compartments_objects_names
    warning('CellOrganizer:instance2VCML', 'Assuming there is only one cell')
    EC_object_name = '';
    cell_object_name = '';
    nucleus_object_name = '';
    % Find EC, cell, and nucleus objects
    for i = 1:n_compartments_objects
        compartment_object_name = all_compartments_objects_names{i};
        compartment_object = all_compartments_objects(compartment_object_name);
        compartment_object_compartment_name = compartment_object.compartment_name;
        if strcmp(compartment_object_compartment_name, EC_translated)
            EC_object_name = compartment_object_name;
        elseif strcmp(compartment_object_compartment_name, cell_translated)
            cell_object_name = compartment_object_name;
        elseif strcmp(compartment_object_compartment_name, nucleus_translated)
            nucleus_object_name = compartment_object_name;
        end
    end
    framework_compartments_object_indices = containers.Map('KeyType', 'char', 'ValueType', 'double');
    framework_compartments_object_indices(EC_translated) = all_compartments_objects_indices(EC_object_name);
    framework_compartments_object_indices(cell_translated) = all_compartments_objects_indices(cell_object_name);
    framework_compartments_object_indices(nucleus_translated) = all_compartments_objects_indices(nucleus_object_name);
    
    object_adjacent_pairs(end+1, 1:2) = [all_compartments_objects(EC_object_name).object_index, all_compartments_objects(cell_object_name).object_index];
    object_adjacent_pairs(end+1, 1:2) = [all_compartments_objects(cell_object_name).object_index, all_compartments_objects(nucleus_object_name).object_index];
    
    for i = 1:n_compartments_objects
        compartment_object_name = all_compartments_objects_names{i};
        compartment_object = all_compartments_objects(compartment_object_name);
        compartment_object_compartment_name = compartment_object.compartment_name;
        compartment_object_volume = compartment_object.Volume;
        compartment_object_index = compartment_object.object_index;
        compartment_object_outside_outside = '';
        compartment_object_outside = network_info.compartments(compartment_object_compartment_name).outside;
        
        if network_info.compartments.isKey(compartment_object_outside)
            compartment_object_outside_outside = network_info.compartments(compartment_object_outside).outside;
        end
        
        if ~is_framework_compartment_name_function(compartment_object_compartment_name)
            if is_framework_compartment_name_function(compartment_object_outside_outside)
                object_adjacent_pairs(end+1, 1:2) = [framework_compartments_object_indices(compartment_object_outside_outside), compartment_object_index];
            end
        end
    end
    
    object_adjacent_pairs = sort(object_adjacent_pairs, 2);
    object_adjacent_pairs = unique(object_adjacent_pairs, 'rows');
end



warning('CellOrganizer:instance2VCML', 'Some inconsistencies! Voxel not 6-connected then determined to be 6-adjacent')


% If was_given_vcml_file is true, instance2VCML will assume VCMLfile is complete and only attempt to replace geometry
was_given_net_file = ~isempty(NETfile);
was_given_vcml_file = ~isempty(VCMLfile);
was_given_biochemistry_file = was_given_net_file || was_given_vcml_file;
if was_given_vcml_file && was_given_net_file
    error('Giving a NET file and a VCML file together not supported');
end

if was_given_vcml_file
    %Create initial Node Object
    docNode = com.mathworks.xml.XMLUtils.createDocument('vcml');
    docNode = xmlread(VCMLfile);
    docRootNode = docNode.getDocumentElement();
    docRootNode = XMLremoveWhitespaceNodes(docRootNode);
    
    if VCMLNumSimulations ~= 1
        error('Not implemented');
    end
    
    warning('CellOrganizer:instance2VCML', 'Make values compatible with unit system or make assumptions!');
else
    %Create initial Node Object
    docNode = com.mathworks.xml.XMLUtils.createDocument('vcml');
    %Create Root node and define the namespace for VCML
    docRootNode = docNode.getDocumentElement();

    docRootNode.setAttribute('xmlns','http://sourceforge.net/projects/vcell/vcml');
    docRootNode.setAttribute('Version', 'Rel_Version_7.1.0_build_3');
end

BioModelNode = createOrGetChild(docRootNode, 'BioModel', [], struct('Name', 'BioModel1'));
BioModelNodeName = char(BioModelNode.getAttribute('Name'));

ModelNode = createOrGetChild(BioModelNode, 'Model', [], struct('Name', 'Model1'));



% Following based on a VCML file saved using Virtual Cell



% Data required to generate VCML

all_compartment_VCML_data = struct();
for i = 1:n_compartments
    name = all_compartment_names{i};
    object_volume = named_imgs_volumes(name);
    all_compartment_VCML_data(i).name = name;
    all_compartment_VCML_data(i).object_volume = object_volume;
    all_compartment_VCML_data(i).FeatureName = name;
    all_compartment_VCML_data(i).DiagramNodeName = all_compartment_VCML_data(i).FeatureName;
    all_compartment_VCML_data(i).PixelClassNodeName = name;
    all_compartment_VCML_data(i).SubVolumeNodeName = name;
end


all_object_VCML_data = struct();
for i = 1:n_compartments_objects
    name = all_compartments_objects_names{i};
    object = all_compartments_objects(name);
    object_volume = object.Volume;
    all_object_VCML_data(i).name = name;
    all_object_VCML_data(i).object = object;
    all_object_VCML_data(i).object_volume = object_volume;
    all_object_VCML_data(i).SpatialObjectName = object.SpatialObjectName;
    all_object_VCML_data(i).VolumeRegionNodeName = name;
    all_object_VCML_data(i).SubVolumeNodeName = all_compartments_objects_names_to_compartment_names(name);
end

all_object_membrane_VCML_data = struct();
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
    all_object_membrane_VCML_data(k).ij = ij;
    all_object_membrane_VCML_data(k).i = i;
    all_object_membrane_VCML_data(k).j = j;
    all_object_membrane_VCML_data(k).name = name;
    all_object_membrane_VCML_data(k).name_i = name_i;
    all_object_membrane_VCML_data(k).name_j = name_j;
    all_object_membrane_VCML_data(k).object_i = object_i;
    all_object_membrane_VCML_data(k).object_j = object_j;
    all_object_membrane_VCML_data(k).object_i_volume = object_i_volume;
    all_object_membrane_VCML_data(k).object_j_volume = object_j_volume;
    all_object_membrane_VCML_data(k).object_i_area = object_i_area;
    all_object_membrane_VCML_data(k).object_j_area = object_j_area;
    all_object_membrane_VCML_data(k).object_area = membrane_area;
    all_object_membrane_VCML_data(k).MembraneRegionNodeName = ['membrane_', name];
    all_object_membrane_VCML_data(k).SpatialObjectName = ['sobj_', name];
    if object_i_volume > object_j_volume
        all_object_membrane_VCML_data(k).subVolumeOutside = all_compartments_objects_names_to_compartment_names(name_i);
        all_object_membrane_VCML_data(k).subVolumeInside = all_compartments_objects_names_to_compartment_names(name_j);
        all_object_membrane_VCML_data(k).regionIdOutside = i-1;
        all_object_membrane_VCML_data(k).regionIdInside = j-1;
    else
        all_object_membrane_VCML_data(k).subVolumeOutside = all_compartments_objects_names_to_compartment_names(name_j);
        all_object_membrane_VCML_data(k).subVolumeInside = all_compartments_objects_names_to_compartment_names(name_i);
        all_object_membrane_VCML_data(k).regionIdOutside = j-1;
        all_object_membrane_VCML_data(k).regionIdInside = i-1;
    end
end


all_membrane_VCML_data = struct();
k = 1;
ball_se = connectivity_se;
for ij = adjacent_pairs'
    i = ij(1);
    j = ij(2);
    start_obj = all_compartment_VCML_data(i);
    link_obj = all_compartment_VCML_data(j);
    
    if start_obj.object_volume > link_obj.object_volume
        outside_obj = start_obj;
        inside_obj = link_obj;
    else
        outside_obj = link_obj;
        inside_obj = start_obj;
    end
    
    % Sort to make independent of order
    [pair_string_translated, pair_string_cell, pair_string] = get_compartment_pair_string(start_obj.name, link_obj.name, name_map);
    
    % Estimate surface area from dilation overlap
    start_image_dilated = single(imdilate(all_compartments_image == i, ball_se));
    link_image_dilated = single(imdilate(all_compartments_image == j, ball_se));
    area = start_image_dilated & link_image_dilated;
    area = sum(area(:)) / 2;
    area = area * prod(resolution);
    
    all_membrane_VCML_data(k).i = i;
    all_membrane_VCML_data(k).j = j;
    all_membrane_VCML_data(k).outside_object = outside_obj;
    all_membrane_VCML_data(k).inside_object = inside_obj;
    all_membrane_VCML_data(k).pair_string = pair_string;
    all_membrane_VCML_data(k).name = pair_string_translated;
    all_membrane_VCML_data(k).area = area;
    all_membrane_VCML_data(k).MembraneNodeName = all_membrane_VCML_data(k).name;
    all_membrane_VCML_data(k).MembraneNodeVoltageName = ['Voltage_', all_membrane_VCML_data(k).name];
    all_membrane_VCML_data(k).DiagramNodeName = all_membrane_VCML_data(k).name;
    all_membrane_VCML_data(k).SurfaceClassNodeName = [all_membrane_VCML_data(k).MembraneNodeName, '_membrane'];
    all_membrane_VCML_data(k).MembraneSubDomainNodeName = all_membrane_VCML_data(k).SurfaceClassNodeName;
    all_membrane_VCML_data(k).SubVolume1Ref = pair_string_cell{1};
    all_membrane_VCML_data(k).SubVolume2Ref = pair_string_cell{2};
    all_membrane_VCML_data(k).OutsideCompartment = pair_string_cell{1};
    all_membrane_VCML_data(k).InsideCompartment = pair_string_cell{2};
    
    k = k + 1;
end
all_membrane_names = {all_membrane_VCML_data.name};
function result2 = is_membrane_function(given_compartments)
    if ischar(given_compartments)
        given_compartments = {given_compartments};
    end
    result2 = cellfun(@(x)any(strcmp(x, all_membrane_names)), given_compartments);
end


% Print for use in BNGL files
fprintf('\n');
fprintf('\n');
fprintf('Membrane areas:\n');
for k = 1:length(all_membrane_VCML_data)
    object = all_membrane_VCML_data(k);
    name = object.name;
    area = object.area;
    fprintf('    %s: %.6f %s\n', name, area, vcml_area_unit);
end
fprintf('\n');
fprintf('\n');


% Add species names
species_digits = 0;
species_name_format = '';
function result2 = species_index_to_name_function(given_species_index)
    result2 = sprintf(species_name_format, given_species_index);
end
function updateSpeciesNames()
    maximum_number_species = length(network_info.species);
    if VCMLAddTranslocationIntermediates
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

parameters_names_to_expressions_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
function result2 = expression_evaluation_function(given_expression, given_return_double)
    if nargin < 2
        return_double = true;
    end
    result2 = mathEval(given_expression, parameters_names_to_expressions_map, struct('return_double', given_return_double));
end

parameter_objects = struct('index', {}, 'name', {}, 'value', {}, 'is_numeric', {}, 'comment', {}, 'units', {});
if was_given_net_file
    updateSpeciesNames();

    % Add or replace volume and surface area parameters
    network_info_parameters_keys = network_info.parameters.keys;
    warning('CellOrganizer:instance2VCML', 'Assuming all unknown parameters have no units');
    for i = 1:length(network_info.parameters)
        parameter_name = network_info_parameters_keys{i};
        parameter = network_info.parameters(parameter_name);
        parameter.units = '1';
        network_info.parameters(parameter_name) = parameter;
    end
    parameter_struct = struct('index', nan, 'name', 'eff_width', 'value_expression', num2str(vcml_effective_width), 'value_expression_parameters', {{}}, 'comment', '', 'units', {vcml_length_unit});
    network_info.parameters(parameter_struct.name) = parameter_struct;
    for i = 1:length(all_compartment_VCML_data)
        object = all_compartment_VCML_data(i);
        name = object.name;
        object_volume = object.object_volume;
        units = vcml_volume_unit;
        parameter_struct = struct('index', nan, 'name', [compartment_volume_prefix,name], 'value_expression', num2str(object_volume), 'value_expression_parameters', {{}}, 'comment', '', 'units', {units});
        network_info.parameters(parameter_struct.name) = parameter_struct;
    end
    for k = 1:length(all_membrane_VCML_data)
        object = all_membrane_VCML_data(k);
        name = object.name;
        area = object.area;
        area_parameter_name = [compartment_area_prefix,name];
        parameter_struct = struct('index', nan, 'name', area_parameter_name, 'value_expression', num2str(area), 'value_expression_parameters', {{}}, 'comment', '', 'units', {vcml_area_unit});
        network_info.parameters(parameter_struct.name) = parameter_struct;
        parameter_struct = struct('index', nan, 'name', [compartment_volume_prefix,name], 'value_expression', [area_parameter_name,'*eff_width'], 'value_expression_parameters', {{area_parameter_name, 'eff_width'}}, 'comment', '', 'units', {vcml_volume_unit});
        network_info.parameters(parameter_struct.name) = parameter_struct;
    end

    % Use a default diffusion coefficient value
    for i = 1:length(network_info.species)
        single_species = network_info.species(i);
        if isnan(single_species.diffusion_coefficient)
            single_species.diffusion_coefficient = vcml_diffusion_coefficient;
        end
        network_info.species(i) = single_species;
    end

    % Set diffusion coefficient units
    network_info_parameters_keys = network_info.parameters.keys;
    warning('CellOrganizer:instance2VCML', 'Assuming dc_* parameters are diffusion coefficients');
    for i = 1:length(network_info.species)
        single_species = network_info.species(i);
        name = single_species.name;
        diffusion_coefficient = single_species.diffusion_coefficient;
        parameter_struct = struct('index', nan, 'name', [compartment_diffusion_coefficient_prefix,name], 'value_expression', num2str(diffusion_coefficient), 'value_expression_parameters', {{}}, 'comment', '', 'units', {vcml_diffusion_coefficient_unit});
        network_info.parameters(parameter_struct.name) = parameter_struct;
        
        single_species.diffusion_coefficient = parameter_struct.name;
        network_info.species(i) = single_species;
    end

    % Set reaction rate constant units
    network_info_parameters_keys = network_info.parameters.keys;
    warning('CellOrganizer:instance2VCML', 'Assuming k_*, kp_*, km_* parameters are reaction rate constants');
    for i = 1:length(network_info.parameters)
        parameter_name = network_info_parameters_keys{i};
        parameter = network_info.parameters(parameter_name);
        if ~((length(parameter_name) > 2 && strcmp(parameter_name(1:2), 'k_')) || (length(parameter_name) > 3 && any(strcmp(parameter_name(1:3), {'kp_', 'km_'}))))
            continue;
        end
        
        % Find reaction associated with constant
        for j = 1:length(network_info.reactions)
            reaction = network_info.reactions(j);
            if ~any(strcmp(parameter_name, reaction.rate_constant_parameters))
                continue;
            end
            % "In general, any nth order volume-independent macroscopic reaction rate constant should have the units molar^{1-n} per second. However, since BNG treats the rate as microscopic, the modeler must convert the macroscopic rate constants to microscopic ones, as follows."
            reaction_n_reactants = length(reaction.reactant_indices);
            if reaction_n_reactants < 1
                error('reaction_n_reactants = %0i', reaction_n_reactants);
            end
            reaction_rate_constant_units = [strjoin([repmat({[vcml_concentration_unit,'-1']},1,reaction_n_reactants-1),{[vcml_time_unit,'-1']}],'.')];
            parameter.units = reaction_rate_constant_units;
            break;
        end
        network_info.parameters(parameter_name) = parameter;
    end

    % Set units for concentration parameters
    network_info_parameters_keys = network_info.parameters.keys';
    network_info_parameters_length = length(network_info.parameters);
    for i = 1:length(network_info.species)
        single_species = network_info.species(i);
        single_species_concentration = single_species.concentration;
        if isStringNumeric(single_species_concentration)
            continue;
        end
        single_species_concentration_identifiers = regexp(single_species_concentration, network_info.regexp_patterns.bng_patterns.string_pattern, 'match');
        single_species_concentration_expressions = regexp(single_species_concentration, network_info.regexp_patterns.bng_patterns.net_expression, 'match');
        single_species_concentration_is_identifier = length(single_species_concentration_identifiers) == 1 && strcmp(single_species_concentration, single_species_concentration_identifiers{1});
        if single_species_concentration_is_identifier
            parameter_name = single_species_concentration;
            parameter = network_info.parameters(parameter_name);
            parameter.units = net_concentration_unit;
            network_info.parameters(parameter_name) = parameter;
        else
            parameter_name = sprintf('Species'); 
            parameter = struct();
            error('Unfinished');
        end
    end
    network_info_parameters_keys = network_info.parameters.keys';
    network_info_parameters_length = length(network_info.parameters)

    % Set initial concentration and count units
    network_info_parameters_keys = network_info.parameters.keys;
    warning('CellOrganizer:instance2VCML', 'Assuming *_init_count parameters have units ''%s''', net_count_unit);
    warning('CellOrganizer:instance2VCML', 'Assuming *_init_conc parameters have units ''%s''', net_concentration_unit);
    warning('CellOrganizer:instance2VCML', 'Assuming N_A parameter is Avogadro constant, has value ''%s'', and has units ''%s''', avogadro_constant_value_expression, avogadro_constant_units);
    for i = 1:length(network_info.parameters)
        parameter_name = network_info_parameters_keys{i};
        parameter = network_info.parameters(parameter_name);
        parameter_name_is_count = strcmp(regexp(parameter_name, ['^',network_info.regexp_patterns.bng_patterns.string_pattern,'_init_count$'], 'match'), parameter_name);
        parameter_name_is_conc = strcmp(regexp(parameter_name, ['^',network_info.regexp_patterns.bng_patterns.string_pattern,'_init_conc$'], 'match'), parameter_name);
        parameter_name_is_avogadro_constant = strcmp(regexp(parameter_name, ['^N_A$'], 'match'), parameter_name);
        parameter_name_is_area = strcmp(regexp(parameter_name, ['^',network_info.regexp_patterns.bng_patterns.string_pattern,'_vol_sphere_area$'], 'match'), parameter_name);
        if parameter_name_is_count
            parameter.units = net_count_unit;
        elseif parameter_name_is_conc
            parameter.units = net_concentration_unit;
        elseif parameter_name_is_avogadro_constant
            parameter.value_expression = avogadro_constant_value_expression;
            parameter.units = avogadro_constant_units;
        elseif parameter_name_is_area
            parameter.units = net_area_unit;
        end
        network_info.parameters(parameter_name) = parameter;
    end

    % Determine if parameters have numeric values
    network_info_parameters_keys = network_info.parameters.keys;
    for i = 1:length(network_info.parameters)
        parameter_name = network_info_parameters_keys{i};
        parameter = network_info.parameters(parameter_name);
        parameter.is_numeric = isStringNumeric(parameter.value_expression);
        network_info.parameters(parameter_name) = parameter;
    end

    network_info_parameters_keys = network_info.parameters.keys';
    network_info_parameters_values = cell2mat(network_info.parameters.values');
    network_info_parameters_values_index = [network_info_parameters_values.index]';
    network_info_parameters_values_name = {network_info_parameters_values.name}';
    network_info_parameters_values_value_expression = {network_info_parameters_values.value_expression}';
    network_info_parameters_values_value_expression_is_numeric = cellfun(@(x)isStringNumeric(x), network_info_parameters_values_value_expression);
    network_info_parameters_values_comment = {network_info_parameters_values.comment}';
    network_info_parameters_values_units = {network_info_parameters_values.units}';

    parameters_names_to_expressions_map = containers.Map(network_info_parameters_values_name, network_info_parameters_values_value_expression);
    parameter_objects = struct('index', network_info_parameters_values_index, 'name', network_info_parameters_values_name, 'value', network_info_parameters_values_value_expression, 'is_numeric', num2cell(network_info_parameters_values_value_expression_is_numeric), 'comment', network_info_parameters_values_comment, 'units', network_info_parameters_values_units);


    % Print for use in BNGL files
    fprintf('\n');
    fprintf('\n');
    fprintf('Parameter values:\n');
    parameter_names_to_print = {};
    parameter_names_to_print{end+1} = 'eff_width';
    for i = 1:length(parameter_names_to_print)
        name = parameter_names_to_print{i};
        object = network_info.parameters(name);
        value_expression = object.value_expression;
        units = object.units;
        value = expression_evaluation_function(value_expression, true);
        fprintf('    %s: %s = %.6f %s\n', name, value_expression, value, units);
    end
    fprintf('\n');
    fprintf('\n');


end



% Set up VCML structure


DiagramNodesMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
DiagramNodesIndices = containers.Map('KeyType', 'char', 'ValueType', 'double');

warning('CellOrganizer:instance2VCML', 'Ignoring CSGdata.primitiveOnly.');
warning('CellOrganizer:instance2VCML', 'Ordinals not used here, connectivity inferred using adjacency in image.');


warning('CellOrganizer:instance2VCML', 'consensus_compartment: Check if this follows rules in Sekar and Faeder 2012');
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
warning('CellOrganizer:instance2VCML', 'Assuming NET initial values are extensive as recommended for BNG');
species_initial_concentrations = cell(length(network_info.species), 1);
for i = 1:length(network_info.species)
    single_species = network_info.species(i);
    concentration = single_species.concentration;
    compartment_size_expression = network_info.compartments(single_species.chosen_compartment).size_expression;
    
    if isStringNumeric(concentration)
        concentration = str2double(concentration);
    end
    
    if isStringNumeric(compartment_size_expression)
        compartment_size_expression = str2double(compartment_size_expression);
    end
    
    % Don't check if this is a membrane or a volume
    % "Molecules in surfaces are assumed to be restricted to a small volume enveloping the surface (i.e., the surface volume is equal to the surface area multiplied by a surface thickness) provided by the modeler."
    
    if isnumeric(concentration)
        concentration = num2str(concentration);
    end
    if isnumeric(compartment_size_expression)
        compartment_size_expression = num2str(compartment_size_expression);
    end
    % Count to mol
    concentration_conversion_factor = ['1/',num2str(avogadro_constant),''];
    % mol to mol/L
    concentration_conversion_factor = [concentration_conversion_factor,'*1/(',compartment_size_expression,'*',num2str(unit_convert(net_volume_unit,'L')),')'];
    %  mol/L to vcml_concentration_unit
    concentration_conversion_factor = [concentration_conversion_factor,'*',num2str(unit_convert('M',vcml_concentration_unit))];
    concentration = [concentration,'*',concentration_conversion_factor];
    
    species_initial_concentrations{i} = concentration;
end
[network_info.species.concentration] = species_initial_concentrations{:};


% Add reaction names
reactions_digits = 0;
reactions_name_format = '';
function result2 = getReactionName(given_reaction)
    result2 = sprintf(reactions_name_format, given_reaction.index);
end
function result2 = getExtendedReactionName(given_reaction)
    result2 = strjoin({'reactants', strjoin(cellfun(@(x)num2str(x), num2cell(given_reaction.reactant_indices), 'UniformOutput', false), '_'), 'products', strjoin(cellfun(@(x)num2str(x), num2cell(given_reaction.product_indices), 'UniformOutput', false), '_')}, '_');
end
function updateReactionsNames()
    maximum_number_reactions = length(network_info.reactions);
    if VCMLAddTranslocationIntermediates
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
warning('CellOrganizer:instance2VCML', 'network_info.reactions.compartment values not known to be correct for BNG and/or VCell');
function result2 = getReactionConsensusCompartment(given_reaction)
    given_reaction_species_indices = [given_reaction.reactant_indices, given_reaction.product_indices];
    given_reaction_species = arrayfun(@(x)network_info.species(x), given_reaction_species_indices);
    given_reaction_species_compartments = cell(length(given_reaction_species), 1);
    for j = 1:length(given_reaction_species)
        given_reaction_species_compartments{j} = given_reaction_species(j).compartments;
    end
    given_reaction_species_compartments = cat(2, given_reaction_species_compartments{:});
    given_reaction_species_compartments = unique(given_reaction_species_compartments);
    result2 = getConsensusCompartment(given_reaction_species_compartments);
end
reactions_compartments = cell(length(network_info.reactions), 1);
for i = 1:length(network_info.reactions)
    reaction = network_info.reactions(i);
    reactions_compartments{i} = getReactionConsensusCompartment(reaction);
end
[network_info.reactions.compartment] = reactions_compartments{:};

% Decide on the reaction rate and units for each reaction
warning('CellOrganizer:instance2VCML', 'network_info.reactions.rate_constant_units values not known to be correct for BNG and/or VCell');
warning('CellOrganizer:instance2VCML', 'network_info.reactions.rate_constant values are not processed as microscopic');
reactions_rates = cell(length(network_info.reactions), 1);
reactions_rate_constants = cell(length(network_info.reactions), 1);
reactions_rates_units = cell(length(network_info.reactions), 1);
reactions_rates_parameters = cell(length(network_info.reactions), 1);
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
    % But NET files do not appear to include these corrections except in comments
    if reaction_n_reactants > 1
        reaction_rate_constant = [ ...
            reaction_rate_constant,'*1/(', ...
                strjoin(repmat({[compartment_volume_prefix,reaction_compartment]}, ...
                [reaction_n_reactants-1,1]),'*'), ...
                ')'];
        reaction_rate_constant_units = [strjoin(repmat({[vcml_concentration_unit,'-1']},1,reaction_n_reactants-1), '.')];
    elseif reaction_n_reactants == 1
        reaction_rate_constant = reaction_rate_constant;
    else
        error('reaction_n_reactants = %0i', reaction_n_reactants);
    end
    reaction_rate_constant_units = [vcml_time_unit,'-1'];
    reaction_rate = reaction_rate_constant;
    reaction_rate_units = reaction_rate_constant_units;
    if reaction_n_reactants > 0
        for j = reaction.reactant_indices
            reaction_rate = [reaction_rate,'*',network_info.species(j).name];
            reaction_rate_units = [reaction_rate_units,'.',vcml_concentration_unit];
        end
    end
    
    reactions_rates{i} = reaction_rate;
    reactions_rate_constants{i} = reaction_rate_constant;
    reactions_rates_units{i} = reaction_rate_units;
    reactions_rates_parameters{i} = regexp(reaction_rate, network_info.regexp_patterns.rate_constant_parameters_pattern, 'match');
end
[network_info.reactions.rate] = reactions_rates{:};
[network_info.reactions.rate_units] = reactions_rates_units{:};
[network_info.reactions.rate_parameters] = reactions_rates_parameters{:};

% Split BNGL translocations because VCell does not accept reactions between non-adjacent compartments
reactions_to_add = network_info.reactions(1:0);
reactions_to_remove = [];
warning('CellOrganizer:instance2VCML', 'Splitting BNGL translocations with assumptions about diffusion coefficients and reaction rate, see code');
warning('CellOrganizer:instance2VCML', 'Permitting BNGL translocations that are invalid according to Sekar and Faeder 2012, "Rule-Based Modeling of Signal Transduction: A Primer" but used in the published example "journal.pcbi.1004611.s003.bngl" with assumptions about diffusion coefficients and reaction rate, see code');
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
    if length(reaction_species_compartments) < 2
        % Do nothing
    elseif length(reaction_species_compartments) == 2
        if is_network_info_compartment_pair_adjacent(reaction_species_compartments{1}, reaction_species_compartments{2}, 1)
            % Do nothing
        else
            if VCMLAddTranslocationIntermediates
                % Find path of compartments from reaction_species_compartments{1} to reaction_species_compartments{2}
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
                    previous_intermediate_reaction2_rate_parameters = previous_intermediate_reaction2.rate_parameters;
                    intermediate_reaction2_rate = previous_intermediate_reaction2_rate;
                    intermediate_reaction2_rate_parameters = previous_intermediate_reaction2_rate_parameters;
                    % Replace reactants in reaction rate of previous_intermediate_reaction2
                    if isnumeric(intermediate_reaction2_rate)
                        intermediate_reaction2_rate = num2str(intermediate_reaction2_rate);
                    end
                    if length(previous_intermediate_reaction2.reactant_indices) ~= 1
                        intermediate_reaction2_rate = ['(',intermediate_reaction2_rate,')^(1/',num2str(length(previous_intermediate_reaction2.reactant_indices)),')'];
                    end
                    previous_intermediate_reaction2_product = network_info.species(previous_intermediate_reaction2.product_indices(1));
                    for k = 1:length(previous_intermediate_reaction2.reactant_indices)
                        previous_intermediate_reaction2_reactant = network_info.species(previous_intermediate_reaction2.reactant_indices(k));
                        intermediate_reaction2_rate = strrep(intermediate_reaction2_rate, previous_intermediate_reaction2_reactant.name, previous_intermediate_reaction2_product.name);
                        
                        intermediate_reaction2_rate_parameters_matching = strcmp(intermediate_reaction2_rate_parameters, previous_intermediate_reaction2_reactant.name);
                        intermediate_reaction2_rate_parameters{intermediate_reaction2_rate_parameters_matching} = previous_intermediate_reaction2_product.name;
                    end
                    intermediate_reaction2_rate_parameters = unique(intermediate_reaction2_rate_parameters);
                    
                    intermediate_reaction2.rate = intermediate_reaction2_rate;
                    intermediate_reaction2.rate_parameters = intermediate_reaction2_rate_parameters;
                    reactions_to_add(end+1) = intermediate_reaction2;
                    
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


if ~was_given_vcml_file
    % Was commented because not known to be necessary. `ModelParameters` does not appear in a VCML file built entirely in VCell itself. Guessing it is generated by SBML import.
    ModelParametersNode = docNode.createElement('ModelParameters');
    ModelNode.appendChild(ModelParametersNode);
    for i = 1:length(parameter_objects)
        object = parameter_objects(i);
        ParameterNode = docNode.createElement('Parameter');
        ModelParametersNode.appendChild(ParameterNode);
        ParameterNode.setAttribute('Name',object.name);
        ParameterNode.setAttribute('Role','user defined');
        ParameterNode.setAttribute('Unit', object.units);
        ParameterNode.setTextContent(object.value);
    end

    for i = 1:length(network_info.species)
        single_species = network_info.species(i);
        
        CompoundNode = docNode.createElement('Compound');
        CompoundNode.setAttribute('Name',single_species.name);
        ModelNode.appendChild(CompoundNode);
    end
end

for i = 1:length(all_compartment_VCML_data)
    object = all_compartment_VCML_data(i);
    name = object.name;
    
    createOrGetChild(ModelNode, 'Feature', struct('Name', object.FeatureName));
end

for k = 1:length(all_membrane_VCML_data)
    object = all_membrane_VCML_data(k);
    i = object.i;
    j = object.j;
    area = object.area;
    
    start_obj = all_compartment_VCML_data(i);
    link_obj = all_compartment_VCML_data(j);
    
    MembraneNode = createOrGetChild(ModelNode, 'Membrane', struct('Name', object.MembraneNodeName));
    if isempty(getAttributeValue(MembraneNode, 'MembraneVoltage', 'str'))
        MembraneNode.setAttribute('MembraneVoltage',object.MembraneNodeVoltageName);
    end
    ModelNode.appendChild(MembraneNode);
end

if ~was_given_vcml_file
    for i = 1:length(network_info.species)
        single_species = network_info.species(i);
        
        LocalizedCompoundNode = docNode.createElement('LocalizedCompound');
        LocalizedCompoundNode.setAttribute('Name',single_species.name);
        LocalizedCompoundNode.setAttribute('SbmlName',single_species.extended_name);
        LocalizedCompoundNode.setAttribute('CompoundRef',single_species.name);
        LocalizedCompoundNode.setAttribute('Structure',single_species.chosen_compartment);
        LocalizedCompoundNode.setAttribute('OverrideName',mat2str(true));
        ModelNode.appendChild(LocalizedCompoundNode);
        
        AnnotationNode = docNode.createElement('Annotation');
        AnnotationNode.setTextContent(single_species.comment);
        LocalizedCompoundNode.appendChild(AnnotationNode);
    end



    warning('CellOrganizer:instance2VCML', 'Rate constant units converted improperly');
    species_rates = repmat({''}, length(network_info.species), 1);
    reaction_mass_action_rates = repmat({''}, length(network_info.species), 1);
    reaction_microscopic_irreversible_rates = repmat({''}, length(network_info.species), 1);
    for i = 1:length(network_info.reactions)
        reaction = network_info.reactions(i);
        
        SimpleReactionNode = docNode.createElement('SimpleReaction');
        SimpleReactionNode.setAttribute('Name',reaction.name);
        SimpleReactionNode.setAttribute('Structure',reaction.compartment);
        % TODO: Check this
        SimpleReactionNode.setAttribute('Reversible',mat2str(false));
        SimpleReactionNode.setAttribute('FluxOption','MolecularOnly');
        ModelNode.appendChild(SimpleReactionNode);
        
        % Reactant nodes and stoichiometries
        reaction_reactant_indices_unique = unique(reaction.reactant_indices);
        reactant_stoichiometries = arrayfun(@(x)sum(reaction.reactant_indices == x), reaction_reactant_indices_unique);
        reactant_first_indices = arrayfun(@(x)find(reaction.reactant_indices == x, 1, 'first'), reaction_reactant_indices_unique);
        for j = 1:length(reaction_reactant_indices_unique)
            reactant_index = reaction.reactant_indices(reactant_first_indices(j));
            ReactantNode = docNode.createElement('Reactant');
            ReactantNode.setAttribute('LocalizedCompoundRef',network_info.species(reactant_index).name);
            ReactantNode.setAttribute('Stoichiometry',num2str(reactant_stoichiometries(j)));
            SimpleReactionNode.appendChild(ReactantNode);
        end
        
        % Product nodes and stoichiometries
        reaction_product_indices_unique = unique(reaction.product_indices);
        product_stoichiometries = arrayfun(@(x)sum(reaction.product_indices == x), reaction_product_indices_unique);
        product_first_indices = arrayfun(@(x)find(reaction.product_indices == x, 1, 'first'), reaction_product_indices_unique);
        for j = 1:length(reaction_product_indices_unique)
            product_index = reaction.product_indices(product_first_indices(j));
            ProductNode = docNode.createElement('Product');
            ProductNode.setAttribute('LocalizedCompoundRef',network_info.species(product_index).name);
            ProductNode.setAttribute('Stoichiometry',num2str(product_stoichiometries(j)));
            SimpleReactionNode.appendChild(ProductNode);
        end
        
        % Construct species rate expressions
        for j = 1:length(reaction_reactant_indices_unique)
            reactant_index = reaction.reactant_indices(reactant_first_indices(j));
            if isempty(species_rates{reactant_index})
                species_rates{reactant_index} = [species_rates{reactant_index},'-'];
            else
                species_rates{reactant_index} = [species_rates{reactant_index},' - '];
            end
            species_rates{reactant_index} = [species_rates{reactant_index},num2str(reactant_stoichiometries(j)),' * (',reaction.rate,')'];
        end
        for j = 1:length(reaction_product_indices_unique)
            product_index = reaction.product_indices(product_first_indices(j));
            if ~isempty(species_rates{product_index})
                species_rates{product_index} = [species_rates{product_index},' + '];
            end
            species_rates{product_index} = [species_rates{product_index},num2str(product_stoichiometries(j)),' * (',reaction.rate,')'];
        end
        
        % Construct reaction rate expression for mass action kinetics
        reaction_mass_action_rate = '';
        reaction_mass_action_rate = [reaction_mass_action_rate, repmat('(', 1, length(reaction_reactant_indices_unique)+1), 'Kf'];
        for j = 1:length(reaction_reactant_indices_unique)
            reactant_index = reaction.reactant_indices(reactant_first_indices(j));
            reaction_mass_action_rate = [reaction_mass_action_rate, ' * ', network_info.species(reactant_index).name, ')'];
        end
        reaction_mass_action_rate = [reaction_mass_action_rate, ' - '];
        reaction_mass_action_rate = [reaction_mass_action_rate, repmat('(', 1, length(reaction_product_indices_unique)), 'Kr'];
        for j = 1:length(reaction_product_indices_unique)
            product_index = reaction.product_indices(product_first_indices(j));
            reaction_mass_action_rate = [reaction_mass_action_rate, ' * ', network_info.species(product_index).name, ')'];
        end
        reaction_mass_action_rate = [reaction_mass_action_rate, ')'];
        reaction_mass_action_rates{i} = reaction_mass_action_rate;
        
        
        % Construct reaction rate expression for microscopic irreversible kinetics
        reaction_microscopic_irreversible_rate = '';
        reaction_microscopic_irreversible_rate = [reaction_microscopic_irreversible_rate, repmat('(', 1, length(reaction_reactant_indices_unique)), 'Kon'];
        for j = 1:length(reaction_reactant_indices_unique)
            reactant_index = reaction.reactant_indices(reactant_first_indices(j));
            reaction_microscopic_irreversible_rate = [reaction_microscopic_irreversible_rate, repmat([' * ', network_info.species(reactant_index).name], [1, reactant_stoichiometries(j)]), ')'];
        end
        reaction_microscopic_irreversible_rates{i} = reaction_microscopic_irreversible_rate;
        
        
        KineticsNode = docNode.createElement('Kinetics');
        SimpleReactionNode.appendChild(KineticsNode);
        
        
        KineticsNode.setAttribute('KineticsType','MassAction');
        
        ParameterNode = docNode.createElement('Parameter');
        KineticsNode.appendChild(ParameterNode);
        ParameterNode.setAttribute('Role','reaction rate');
        ParameterNode.setAttribute('Unit',reaction.rate_units);
        
        ParameterNode.setAttribute('Name','J');
        ParameterNode.setTextContent(reaction_mass_action_rate);
        
        if is_membrane_function(reaction.compartment)
            reaction_forward_rate_constant_units = [strjoin([repmat({[vcml_concentration_unit,'-1']},1,length(reaction.reactant_indices)-1),{[vcml_time_unit,'-1']}],'.')];
            reaction_reverse_rate_constant_units = [strjoin([repmat({[vcml_membrane_substance_unit,'-1']},1,length(reaction.product_indices)-1),{[vcml_time_unit,'-1']}],'.')];
            
            reaction_rate_units = vcml_membrane_reaction_rate_unit;
        else
            reaction_forward_rate_constant_units = [strjoin([repmat({[vcml_concentration_unit,'-1']},1,length(reaction.reactant_indices)-1),{[vcml_time_unit,'-1']}],'.')];
            reaction_reverse_rate_constant_units = [strjoin([repmat({[vcml_concentration_unit,'-1']},1,length(reaction.product_indices)-1),{[vcml_time_unit,'-1']}],'.')];
            
            reaction_rate_units = vcml_volume_reaction_rate_unit;
        end

        ParameterNode.setAttribute('Unit',reaction_rate_units);
        
        ParameterNode = docNode.createElement('Parameter');
        ParameterNode.setAttribute('Name','Kf');
        ParameterNode.setAttribute('Role','forward rate constant');
        ParameterNode.setAttribute('Unit',reaction_forward_rate_constant_units);
        ParameterNode.setTextContent(reaction.rate_constant);
        KineticsNode.appendChild(ParameterNode);

        ParameterNode = docNode.createElement('Parameter');
        ParameterNode.setAttribute('Name','Kr');
        ParameterNode.setAttribute('Role','reverse rate constant');
        ParameterNode.setAttribute('Unit',reaction_reverse_rate_constant_units);
        ParameterNode.setTextContent(num2str(0));
        KineticsNode.appendChild(ParameterNode);
        
        AnnotationNode = docNode.createElement('Annotation');
        AnnotationNode.setTextContent(reaction.comment);
        SimpleReactionNode.appendChild(AnnotationNode);
    end
    warning('CellOrganizer:instance2VCML', 'ParameterNode with role ''user defined'' attribute Unit not computed, assuming ??? for 1 parameter');
    [network_info.species.rate] = species_rates{:};
    [network_info.reactions.mass_action_rate] = reaction_mass_action_rates{:};
end



for i = 1:length(all_compartment_VCML_data)
    object = all_compartment_VCML_data(i);
    name = object.name;
    
    DiagramNode = createOrGetChild(ModelNode, 'Diagram', struct('Name', object.DiagramNodeName, 'Structure', object.FeatureName));
    
    DiagramNodesMap(object.FeatureName) = DiagramNode;
    DiagramNodesIndices(object.FeatureName) = length(DiagramNodesMap);
end

for k = 1:length(all_membrane_VCML_data)
    object = all_membrane_VCML_data(k);
    i = object.i;
    j = object.j;
    area = object.area;
    
    start_obj = all_compartment_VCML_data(i);
    link_obj = all_compartment_VCML_data(j);
    
    DiagramNode = createOrGetChild(ModelNode, 'Diagram', struct('Name', object.DiagramNodeName, 'Structure', object.MembraneNodeName));
    
    DiagramNodesMap(object.DiagramNodeName) = DiagramNode;
    DiagramNodesIndices(object.DiagramNodeName) = length(DiagramNodesMap);
end

ModelUnitSystemNode_attributes = cell(0, 3);
ModelUnitSystemNode_attributes(end+1, [1, 2]) = {'VolumeSubstanceUnit',vcml_volume_substance_unit};
ModelUnitSystemNode_attributes(end+1, [1, 2]) = {'MembraneSubstanceUnit',vcml_membrane_substance_unit};
ModelUnitSystemNode_attributes(end+1, [1, 2]) = {'LumpedReactionSubstanceUnit','molecules'};
ModelUnitSystemNode_attributes(end+1, [1, 2]) = {'VolumeUnit',vcml_volume_unit};
ModelUnitSystemNode_attributes(end+1, [1, 2]) = {'LengthUnit',vcml_length_unit};
ModelUnitSystemNode_attributes(end+1, [1, 2]) = {'AreaUnit',vcml_area_unit};
ModelUnitSystemNode_attributes(end+1, [1, 2]) = {'TimeUnit',vcml_time_unit};

function result2 = convert_from_vcml_to_given_vcml(given_value, given_unit_name)
    if ~was_given_vcml_file
        result2 = given_value;
        return;
    end
    result2 = [];
    for ii = 1:size(ModelUnitSystemNode_attributes, 1)
        if strcmp(ModelUnitSystemNode_attributes{ii, 1}, given_unit_name)
            result2 = given_value * unit_convert(ModelUnitSystemNode_attributes{ii, 2}, ModelUnitSystemNode_attributes{ii, 3});
            break;
        end
    end
    if isempty(result2)
        error('Could not convert value with given_unit_name "%s"', given_unit_name);
    end
end

function result2 = get_given_vcml_unit(given_unit_name)
    if ~was_given_vcml_file
        result2 = given_unit_name;
        return;
    end
    result2 = [];
    for ii = 1:size(ModelUnitSystemNode_attributes, 1)
        if strcmp(ModelUnitSystemNode_attributes{ii, 1}, given_unit_name)
            result2 = ModelUnitSystemNode_attributes{ii, 3};
            break;
        end
    end
    if isempty(result2)
        error('Could not get unit with given_unit_name "%s"', given_unit_name);
    end
end

should_assert_ModelUnitSystemNode_attributes_equal = false;
if was_given_vcml_file
    warning('CellOrganizer:instance2VCML', 'Not checking ModelUnitSystem');
    ModelUnitSystemNode = ModelNode.getFirstChildByTagName('ModelUnitSystem');
    for i = 1:size(ModelUnitSystemNode_attributes, 1)
        vcml_file_model_unit = javaStringToValue(ModelUnitSystemNode.getAttribute(ModelUnitSystemNode_attributes{i, 1}), 'char');
        ModelUnitSystemNode_attributes{i, 3} = vcml_file_model_unit;
        if should_assert_ModelUnitSystemNode_attributes_equal
            assert(valuesEqual(ModelUnitSystemNode_attributes{i, 2}, ModelUnitSystemNode_attributes{i, 3}));
        end
    end
else
    ModelUnitSystemNode_attributes(:, 3) = ModelUnitSystemNode_attributes(:, 2);
end
if ~was_given_vcml_file
    ModelUnitSystemNode = docNode.createElement('ModelUnitSystem');
    warning('CellOrganizer:instance2VCML', 'Unit conversions unfinished');
    warning('CellOrganizer:instance2VCML', 'Check VolumeSubstanceUnit (''%s'')',vcml_volume_substance_unit);
    for i = 1:size(ModelUnitSystemNode_attributes, 1)
        ModelUnitSystemNode.setAttribute(ModelUnitSystemNode_attributes{i, 1}, ModelUnitSystemNode_attributes{i, 2});
    end
    ModelNode.appendChild(ModelUnitSystemNode);
end


% Add MembraneVoltage constants here so they are in MathDescription Constants but not in ModelParameters Parameters
for j = 1:length(all_membrane_VCML_data)
    object = all_membrane_VCML_data(j);
    parameter = struct('index', nan, 'name', object.MembraneNodeVoltageName, 'value_expression', '0', 'value_expression_parameters', {{}}, 'units', 'V', 'is_numeric', true);
    network_info.parameters(object.MembraneNodeVoltageName) = parameter;
end

% TODO: OutputFunctions?
% TODO: Stochastic spatial simulations?
% TODO: MembraneSubDomain attributes InsideCompartment and OutsideCompartment should be set correctly
% TODO: Needs SimulationContext, etc.? (https://github.com/virtualcell/vcell/blob/master/vcell-core/src/main/java/cbit/vcell/xml/XmlReader.java#L6128)

% Create simulations
for i = 1:VCMLNumSimulations
    SimulationSpecNode = [];
    GeometryNode = [];
    SimulationSpecNodes = BioModelNode.getElementsByTagName('SimulationSpec');
    if was_given_vcml_file
        % Find the first SimulationSpec with Geometry of dimension 3
        for j = 1:SimulationSpecNodes.getLength()
            SimulationSpecNode2 = SimulationSpecNodes.item(j-1);
            GeometryNode2 = getFirstChildByTagName(SimulationSpecNode2, 'Geometry');
            GeometryNode2_dimension = getAttributeValue(GeometryNode2, 'Dimension', 'int');
            if GeometryNode2_dimension == 3
                SimulationSpecNode = SimulationSpecNode2;
                GeometryNode = GeometryNode2;
                break;
            end
        end
    end
    if ~isempty(SimulationSpecNode)
        SimulationSpecNodeName = char(SimulationSpecNode.getAttribute('Name'));
    else
        SimulationSpecNodeName = sprintf('Spatial%i', i);
    end
    if isempty(SimulationSpecNode)
        SimulationSpecNode = docNode.createElement('SimulationSpec');
        SimulationSpecNode.setAttribute('Name',SimulationSpecNodeName);
        SimulationSpecNode.setAttribute('InsufficientIterations','false');
        SimulationSpecNode.setAttribute('InsufficientMaxMolecules','false');
        SimulationSpecNode.setAttribute('CharacteristicSize','1.0');
        BioModelNode.appendChild(SimulationSpecNode);
        
        NetworkConstraintsNodes = BioModelNode.getElementsByTagName('NetworkConstraints');
        if NetworkConstraintsNodes.getLength() > 0
            NetworkConstraintsNode = NetworkConstraintsNodes.item(0);
        else
            NetworkConstraintsNode = docNode.createElement('NetworkConstraints');
            NetworkConstraintsNode.setAttribute('RbmMaxIteration','3');
            NetworkConstraintsNode.setAttribute('RbmMaxMoleculesPerSpecies','10');
            NetworkConstraintsNode.setAttribute('RbmSpeciesLimit','800');
            NetworkConstraintsNode.setAttribute('RbmReactionsLimit','2500');
            SimulationSpecNode.appendChild(NetworkConstraintsNode);
        end
        SimulationSpecNode.setAttribute('Stochastic','false');
        % SimulationSpecNode.setAttribute('Stochastic','true');
        SimulationSpecNode.setAttribute('UseConcentration','true');
        SimulationSpecNode.setAttribute('RuleBased','false');
    end
    
    SimulationSpecNodes = BioModelNode.getElementsByTagName('SimulationSpec');
    for j = 1:SimulationSpecNodes.getLength()
        SimulationSpecNode2 = SimulationSpecNodes.item(j-1);
        SimulationSpecNode2Name = char(SimulationSpecNode2.getAttribute('Name'));
        SimulationSpecNode2Name = [options.prefix, '_', SimulationSpecNode2Name];
        SimulationSpecNode2.setAttribute('Name',SimulationSpecNode2Name);
        SimulationNodes2 = SimulationSpecNode2.getElementsByTagName('Simulation');
        for k = 1:SimulationNodes2.getLength()
            SimulationNode2 = SimulationNodes2.item(k-1);
            SimulationSpecNode2Name = char(SimulationSpecNode2.getAttribute('Name'));
            SimulationNode2Name = char(SimulationNode2.getAttribute('Name'));
            SimulationNode2Name = [SimulationSpecNode2Name, '_', SimulationNode2Name];
            SimulationNode2.setAttribute('Name',SimulationNode2Name);
        end
    end
    SimulationSpecNodeName = char(SimulationSpecNode.getAttribute('Name'));
    SimulationNodeName = [SimulationSpecNodeName, '_Simulation_', num2str(i)];
    
    
    % Geometry
    
    if was_given_vcml_file
        GeometryNode = getFirstChildByTagName(SimulationSpecNode, 'Geometry');
    else
        GeometryNode = docNode.createElement('Geometry');
        GeometryNodeName = sprintf('Geometry%i', i);
        GeometryNode.setAttribute('Name',GeometryNodeName);
        SimulationSpecNode.appendChild(GeometryNode);
    end
    GeometryNode.setAttribute('Dimension','3');
    
    ExtentNode = createOrGetChild(GeometryNode, 'Extent');
    named_img_key = named_imgs.keys();
    named_img_key = named_img_key{1};
    named_img = named_imgs(named_img_key);
    extent_ij = size(named_img) .* resolution;
    extent_ij_given_vcml = convert_from_vcml_to_given_vcml(extent_ij, 'LengthUnit');
    ExtentNode.setAttribute('X',num2str(extent_ij_given_vcml(2)));
    ExtentNode.setAttribute('Y',num2str(extent_ij_given_vcml(1)));
    ExtentNode.setAttribute('Z',num2str(extent_ij_given_vcml(3)));
    
    OriginNode = createOrGetChild(GeometryNode, 'Origin');
    OriginNode.setAttribute('X',num2str(0));
    OriginNode.setAttribute('Y',num2str(0));
    OriginNode.setAttribute('Z',num2str(0));
    
    ImageNode = createOrGetChild(GeometryNode, 'Image');
    ImageNodeName = 'Image1';
    ImageNode.setAttribute('Name',ImageNodeName);
    
    % https://github.com/virtualcell/vcell/blob/master/vcell-core/src/main/java/cbit/vcell/xml/Xmlproducer.java creates a `VCImageCompressed`, which uses InflaterInputStream to decompress pixel data.
    SamplesBytes = all_compartments_image;
    SamplesBytes = permute(SamplesBytes, [2, 1, 3]);
    SamplesBytes = SamplesBytes(:);
    
    % https://undocumentedmatlab.com/blog/savezip-utility
    inputStream = uint8(SamplesBytes);
    outputStream = java.io.ByteArrayOutputStream();
    compressedInputStream = java.util.zip.DeflaterOutputStream(outputStream);
    compressedInputStream.write(inputStream, 0, numel(inputStream));
    clear inputStream;
    compressedInputStream.finish();
    compressedInputStream.close();
    outputStream.close();
    import javax.xml.bind.DatatypeConverter;
    SamplesBytes = reshape(outputStream.toByteArray(), 1, []);
    SamplesText = DatatypeConverter.printHexBinary(SamplesBytes);
    SamplesText = char(SamplesText);
    SamplesTextLength = length(SamplesText);
    
    ImageDataNode = createOrGetChild(ImageNode, 'ImageData');
    ImageDataNode.setAttribute('X',num2str(size(all_compartments_image,2)));
    ImageDataNode.setAttribute('Y',num2str(size(all_compartments_image,1)));
    ImageDataNode.setAttribute('Z',num2str(size(all_compartments_image,3)));
    ImageDataNode.setAttribute('CompressedSize',num2str(SamplesTextLength));
    ImageDataNode.setTextContent(SamplesText);
    clear SamplesText;
    
    removeChildrenByTagName(ImageNode, 'PixelClass');
    for j = 1:length(all_compartment_VCML_data)
        object = all_compartment_VCML_data(j);
        PixelClassNode = docNode.createElement('PixelClass');
        PixelClassNode.setAttribute('Name',object.PixelClassNodeName);
        PixelClassNode.setAttribute('ImagePixelValue',num2str(j));
        ImageNode.appendChild(PixelClassNode);
    end
    
    removeChildrenByTagName(GeometryNode, 'SubVolume');
    for j = 1:length(all_compartment_VCML_data)
        object = all_compartment_VCML_data(j);
        SubVolumeNode = docNode.createElement('SubVolume');
        SubVolumeNode.setAttribute('Name',object.SubVolumeNodeName);
        SubVolumeNode.setAttribute('Handle',num2str(j-1));
        SubVolumeNode.setAttribute('Type','Image');
        SubVolumeNode.setAttribute('ImagePixelValue',num2str(j));
        GeometryNode.appendChild(SubVolumeNode);
    end
    
    removeChildrenByTagName(GeometryNode, 'SurfaceClass');
    for j = 1:length(all_membrane_VCML_data)
        object = all_membrane_VCML_data(j);
        SurfaceClassNode = docNode.createElement('SurfaceClass');
        SurfaceClassNode.setAttribute('Name',object.SurfaceClassNodeName);
        SurfaceClassNode.setAttribute('SubVolume1Ref',object.SubVolume1Ref);
        SurfaceClassNode.setAttribute('SubVolume2Ref',object.SubVolume2Ref);
        GeometryNode.appendChild(SurfaceClassNode);
    end
    
    removeChildrenByTagName(GeometryNode, 'SurfaceDescription');
    SurfaceDescriptionNode = docNode.createElement('SurfaceDescription');
    SurfaceDescriptionNode.setAttribute('NumSamplesX',num2str(size(all_compartments_image,2)));
    SurfaceDescriptionNode.setAttribute('NumSamplesY',num2str(size(all_compartments_image,1)));
    SurfaceDescriptionNode.setAttribute('NumSamplesZ',num2str(size(all_compartments_image,3)));
    SurfaceDescriptionNode.setAttribute('CutoffFrequency','0.3');
    GeometryNode.appendChild(SurfaceDescriptionNode);
    
    for j = 1:n_compartments_objects
        object = all_object_VCML_data(j);
        object_volume = object.object_volume;
        object_volume_given_vcml = convert_from_vcml_to_given_vcml(object_volume, 'VolumeUnit');
        vcml_volume_unit_given_vcml = get_given_vcml_unit('VolumeUnit');
        VolumeRegionNode = docNode.createElement('VolumeRegion');
        VolumeRegionNode.setAttribute('Name',object.VolumeRegionNodeName);
        VolumeRegionNode.setAttribute('RegionID',num2str(j-1));
        VolumeRegionNode.setAttribute('SubVolume',object.SubVolumeNodeName);
        VolumeRegionNode.setAttribute('Size',num2str(object_volume_given_vcml));
        VolumeRegionNode.setAttribute('Unit',vcml_volume_unit_given_vcml);
        SurfaceDescriptionNode.appendChild(VolumeRegionNode);
    end
    
    for j = 1:length(all_object_membrane_VCML_data)
        object = all_object_membrane_VCML_data(j);
        object_area = object.object_area;
        object_area_given_vcml = convert_from_vcml_to_given_vcml(object_area, 'AreaUnit');
        vcml_area_unit_given_vcml = get_given_vcml_unit('AreaUnit');
        ij = object.ij;
        start_obj = all_object_VCML_data(ij(1));
        link_obj = all_object_VCML_data(ij(2));
        
        MembraneRegionNode = docNode.createElement('MembraneRegion');
        MembraneRegionNode.setAttribute('Name',object.MembraneRegionNodeName);
        MembraneRegionNode.setAttribute('VolumeRegion1',start_obj.VolumeRegionNodeName);
        MembraneRegionNode.setAttribute('VolumeRegion2',link_obj.VolumeRegionNodeName);
        MembraneRegionNode.setAttribute('Size',num2str(object_area_given_vcml));
        MembraneRegionNode.setAttribute('Unit',vcml_area_unit_given_vcml);
        SurfaceDescriptionNode.appendChild(MembraneRegionNode);
    end
    
    
    % GeometryContext
    
    removeChildrenByTagName(SimulationSpecNode, 'GeometryContext');
    GeometryContextNode = docNode.createElement('GeometryContext');
    SimulationSpecNode.appendChild(GeometryContextNode);
    
    
    % Assume all other SimulationSpec nodes are non-spatial duplicates and modify their geometric parameters
    for j = 1:SimulationSpecNodes.getLength()
        SimulationSpecNode2 = SimulationSpecNodes.item(j-1);
        GeometryNode2 = getFirstChildByTagName(SimulationSpecNode2, 'Geometry');
        GeometryContextNode2 = getFirstChildByTagName(SimulationSpecNode2, 'GeometryContext');
        GeometryNode2_dimension = getAttributeValue(GeometryNode2, 'Dimension', 'int');
        if GeometryNode2_dimension == 0
            volume_name = getAttributeValue(getFirstChildByTagName(GeometryNode2, 'SubVolume'), 'Name', 'str');
            surface_name = volume_name;
        end
        
        for k = 1:length(all_compartment_VCML_data)
            object = all_compartment_VCML_data(k);
            object_volume = object.object_volume;
            object_volume_given_vcml = convert_from_vcml_to_given_vcml(object_volume, 'VolumeUnit');
            if GeometryNode2_dimension == 3
                volume_name = object.SubVolumeNodeName;
            end
            
            FeatureMappingNode = createOrGetChild(GeometryContextNode2, 'FeatureMapping', struct('Feature', object.FeatureName));
            FeatureMappingNode.setAttribute('Feature',object.FeatureName);
            FeatureMappingNode.setAttribute('GeometryClass',volume_name);
            FeatureMappingNode.setAttribute('SubVolume',volume_name);
            FeatureMappingNode.setAttribute('Size',num2str(object_volume_given_vcml));
            FeatureMappingNode.setAttribute('VolumePerUnitVolume',num2str(1));
            
            BoundariesTypesNode = createOrGetChild(FeatureMappingNode, 'BoundariesTypes');
            BoundariesTypesNode.setAttribute('Xm','Flux');
            BoundariesTypesNode.setAttribute('Xp','Flux');
            BoundariesTypesNode.setAttribute('Ym','Flux');
            BoundariesTypesNode.setAttribute('Yp','Flux');
            BoundariesTypesNode.setAttribute('Zm','Flux');
            BoundariesTypesNode.setAttribute('Zp','Flux');
        end
    
        for k = 1:length(all_membrane_VCML_data)
            object = all_membrane_VCML_data(k);
            object_area = object.area;
            object_area_given_vcml = convert_from_vcml_to_given_vcml(object_area, 'AreaUnit');
            vcml_area_unit_given_vcml = get_given_vcml_unit('AreaUnit');
            if GeometryNode2_dimension == 3
                surface_name = object.SurfaceClassNodeName;
            end
            
            MembraneMappingNode = createOrGetChild(GeometryContextNode2, 'MembraneMapping', struct('Membrane', object.MembraneNodeName));
            MembraneMappingNode.setAttribute('Membrane',object.MembraneNodeName);
            MembraneMappingNode.setAttribute('GeometryClass',surface_name);
            MembraneMappingNode.setAttribute('Size',num2str(object_area_given_vcml));
            MembraneMappingNode.setAttribute('AreaPerUnitArea',num2str(1));
            MembraneMappingNode.setAttribute('CalculateVoltage','false');
            MembraneMappingNode.setAttribute('SpecificCapacitance',num2str(1));
            MembraneMappingNode.setAttribute('InitialVoltage',num2str(0));
        end
    end
    
    
    % ReactionContext
    
    ReactionContextNode = createOrGetChild(SimulationSpecNode, 'ReactionContext');
    
    for j = 1:length(network_info.species)
        object = network_info.species(j);
        object_initial_concentration = object.concentration;
        object_diffusion_coefficient = object.diffusion_coefficient;
        if isnumeric(object_initial_concentration)
            object_initial_concentration = num2str(object_initial_concentration);
        end
        if isnumeric(object_diffusion_coefficient)
            object_diffusion_coefficient = num2str(object_diffusion_coefficient);
        end
        
        LocalizedCompoundSpecNode = docNode.createElement('LocalizedCompoundSpec');
        LocalizedCompoundSpecNode.setAttribute('LocalizedCompoundRef',object.name);
        LocalizedCompoundSpecNode.setAttribute('ForceConstant',mat2str(false));
        LocalizedCompoundSpecNode.setAttribute('WellMixed',mat2str(false));
        LocalizedCompoundSpecNode.setAttribute('ForceContinuous',mat2str(false));
        ReactionContextNode.appendChild(LocalizedCompoundSpecNode);
        
        InitialConcentrationNode = docNode.createElement('InitialConcentration');
        InitialConcentrationNode.setTextContent(object_initial_concentration);
        
        LocalizedCompoundSpecNode.appendChild(InitialConcentrationNode);
        
        DiffusionNode = docNode.createElement('Diffusion');
        DiffusionNode.setTextContent(object_diffusion_coefficient);
        LocalizedCompoundSpecNode.appendChild(DiffusionNode);
    end
    
    for j = 1:length(network_info.reactions)
        reaction = network_info.reactions(j);
        
        ReactionSpecNode = docNode.createElement('ReactionSpec');
        ReactionSpecNode.setAttribute('ReactionStepRef',reaction.name);
        ReactionSpecNode.setAttribute('ReactionMapping','included');
        ReactionContextNode.appendChild(ReactionSpecNode);
    end
    
    
    % MathDescription
    
    MathDescriptionNode = createOrGetChild(SimulationSpecNode, 'MathDescription');
    MathDescriptionNode.setAttribute('Name',[SimulationSpecNodeName, '_generated']);
    
    constant_nodes_to_add = {};
    function_nodes_to_add = {};
    for j = 1:length(parameter_objects)
        parameter_object = parameter_objects(j);
        ParameterNodeName = char(parameter_object.name);
        ParameterNodeValue = parameter_object.value;
        parameter_value_is_numeric = parameter_object.is_numeric;
        if parameter_value_is_numeric
            % Number
            ConstantNode = docNode.createElement('Constant');
            ConstantNode.setAttribute('Name',ParameterNodeName);
            ConstantNode.setTextContent(ParameterNodeValue);
            constant_nodes_to_add{end+1} = ConstantNode;
        else
            % Expression
            FunctionNode = docNode.createElement('Function');
            FunctionNode.setAttribute('Name',ParameterNodeName);
            FunctionNode.setTextContent(ParameterNodeValue);
            function_nodes_to_add{end+1} = FunctionNode;
        end
    end
    
    % Constant before Function
    for j = 1:length(constant_nodes_to_add)
        MathDescriptionNode.appendChild(constant_nodes_to_add{j});
    end
    
    for j = 1:length(function_nodes_to_add)
        MathDescriptionNode.appendChild(function_nodes_to_add{j});
    end
    
    CompartmentSubDomainNodeMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for j = 1:length(all_compartment_VCML_data)
        object = all_compartment_VCML_data(j);
        
        CompartmentSubDomainNode = createOrGetChild(MathDescriptionNode, 'CompartmentSubDomain', struct('Name', object.name));
        CompartmentSubDomainNodeMap(object.name) = CompartmentSubDomainNode;
        
        removeChildrenByTagName(CompartmentSubDomainNode, 'BoundaryType');
        for dim_char = dim_chars
            for sign_char = sign_chars
                BoundaryTypeNode = docNode.createElement('BoundaryType');
                BoundaryTypeNode.setAttribute('Boundary',[dim_char,sign_char]);
                BoundaryTypeNode.setAttribute('Type','Flux');
                CompartmentSubDomainNode.appendChild(BoundaryTypeNode);
            end
        end
    end
    
    
    warning('CellOrganizer:instance2VCML', 'Unfinished')
    
    MembraneSubDomainNodeMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    % Keep track of which MembraneSubDomain nodes were in the input file
    temp = MathDescriptionNode.getElementsByTagName('MembraneSubDomain');
    input_MembraneSubDomainNodes = {};
    for j = 1:temp.getLength()
        input_MembraneSubDomainNodes{end+1} = temp.item(j-1);
    end
    
    for j = 1:length(all_membrane_VCML_data)
        object = all_membrane_VCML_data(j);
        
        MembraneSubDomainNode = createOrGetChild(MathDescriptionNode, 'MembraneSubDomain', struct('InsideCompartment', object.InsideCompartment, 'OutsideCompartment', object.OutsideCompartment));
        MembraneSubDomainNode2 = createOrGetChild(MathDescriptionNode, 'MembraneSubDomain', struct('InsideCompartment', object.OutsideCompartment, 'OutsideCompartment', object.InsideCompartment));
        if ~was_given_biochemistry_file || inCellArray(input_MembraneSubDomainNodes, MembraneSubDomainNode2)
            MathDescriptionNode.removeChild(MembraneSubDomainNode);
            MembraneSubDomainNode = MembraneSubDomainNode2;
        else
            MathDescriptionNode.removeChild(MembraneSubDomainNode2);
        end
        
        MembraneSubDomainNode.setAttribute('Name',object.MembraneSubDomainNodeName);
        MembraneSubDomainNodeMap(object.MembraneSubDomainNodeName) = MembraneSubDomainNode;
        
        removeChildrenByTagName(MembraneSubDomainNode, 'BoundaryType');
        for dim_char = dim_chars
            for sign_char = sign_chars
                BoundaryTypeNode = docNode.createElement('BoundaryType');
                BoundaryTypeNode.setAttribute('Boundary',[dim_char,sign_char]);
                BoundaryTypeNode.setAttribute('Type','Value');
                MembraneSubDomainNode.appendChild(BoundaryTypeNode);
            end
        end
    end
    
    if was_given_biochemistry_file
        % Remove MembraneSubDomain nodes not in generated geometry
        MembraneSubDomainNodeMap_keys = MembraneSubDomainNodeMap.keys();
        for j = 1:length(input_MembraneSubDomainNodes)
            input_MembraneSubDomainNode = input_MembraneSubDomainNodes{j};
            should_remove_input_MembraneSubDomainNode = true;
            for k = 1:length(MembraneSubDomainNodeMap)
                MembraneSubDomainNodeMap_key = MembraneSubDomainNodeMap_keys{k};
                generated_membraneSubDomainNode = MembraneSubDomainNodeMap(MembraneSubDomainNodeMap_key);
                if input_MembraneSubDomainNode == generated_membraneSubDomainNode
                    should_remove_input_MembraneSubDomainNode = false;
                    break;
                end
            end
            if should_remove_input_MembraneSubDomainNode
                MathDescriptionNode.removeChild(input_MembraneSubDomainNode);
            end
        end
    end
    
    
    warning('CellOrganizer:instance2VCML', 'JumpCondition always set to zero');
    warning('CellOrganizer:instance2VCML', 'Unfinished');
    for j = 1:length(network_info.species)
        single_species = network_info.species(j);
        single_species_compartment = single_species.chosen_compartment;
        reaction_rate_constant = reaction.rate_constant;
        
        single_species_initial_concentration = single_species.concentration;
        single_species_diffusion_coefficient = single_species.diffusion_coefficient;
        if isnumeric(single_species_initial_concentration)
            single_species_initial_concentration = num2str(single_species_initial_concentration);
        end
        if isnumeric(single_species_diffusion_coefficient)
            single_species_diffusion_coefficient = num2str(single_species_diffusion_coefficient);
        end
        
        if is_membrane_function(single_species_compartment)
            parent_node = MembraneSubDomainNodeMap(single_species_compartment);
        else
            parent_node = CompartmentSubDomainNodeMap(single_species_compartment);
        end
        
        PdeEquationNode = docNode.createElement('PdeEquation');
        PdeEquationNode.setAttribute('Name',single_species.name);
        PdeEquationNode.setAttribute('SolutionType','Unknown');
        parent_node.appendChild(PdeEquationNode);
        
        if ~is_membrane_function(single_species_compartment)
            for k = 1:length(all_membrane_VCML_data)
                object = all_membrane_VCML_data(k);
                if strcmp(object.InsideCompartment, single_species_compartment) || strcmp(object.OutsideCompartment, single_species_compartment)
                    JumpConditionNode = docNode.createElement('JumpCondition');
                    JumpConditionNode.setAttribute('Name',single_species.name);
                    MembraneSubDomainNodeMap(object.MembraneSubDomainNodeName).appendChild(JumpConditionNode);
                    
                    InFluxNode = docNode.createElement('InFlux');
                    InFluxNode.setTextContent(num2str(0));
                    JumpConditionNode.appendChild(InFluxNode);
                    OutFluxNode = docNode.createElement('OutFlux');
                    OutFluxNode.setTextContent(num2str(0));
                    JumpConditionNode.appendChild(OutFluxNode);
                end
            end
        end
        
        RateNode = docNode.createElement('Rate');
        RateNode.setTextContent(single_species.rate);
        PdeEquationNode.appendChild(RateNode);
        
        DiffusionNode = docNode.createElement('Diffusion');
        DiffusionNode.setTextContent(single_species_diffusion_coefficient);
        PdeEquationNode.appendChild(DiffusionNode);
        
        InitialNode = docNode.createElement('Initial');
        InitialNode.setTextContent(single_species_initial_concentration);
        PdeEquationNode.appendChild(InitialNode);
    end
    
    
    
    % Simulation
    
    removeChildrenByTagName(SimulationSpecNode, 'Simulation');
    
    SimulationNode = docNode.createElement('Simulation');
    SimulationNode.setAttribute('Name',SimulationNodeName);
    SimulationSpecNode.appendChild(SimulationNode);
    
    SolverTaskDescriptionNode = docNode.createElement('SolverTaskDescription');
    SolverTaskDescriptionNode.setAttribute('TaskType','Unsteady');
    SolverTaskDescriptionNode.setAttribute('UseSymbolicJacobian','false');
    SolverTaskDescriptionNode.setAttribute('Solver','Sundials Stiff PDE Solver (Variable Time Step)');
    SimulationNode.appendChild(SolverTaskDescriptionNode);
    
    TimeBoundNode = docNode.createElement('TimeBound');
    TimeBoundNode.setAttribute('StartTime','0.0');
    TimeBoundNode.setAttribute('EndTime',num2str(VCMLEndTime));
    SolverTaskDescriptionNode.appendChild(TimeBoundNode);
    
    TimeStepNode = docNode.createElement('TimeStep');
    TimeStepNode.setAttribute('DefaultTime',num2str(VCMLDefaultTimeStep));
    TimeStepNode.setAttribute('MinTime',num2str(VCMLMinTimeStep));
    TimeStepNode.setAttribute('MaxTime',num2str(VCMLMaxTimeStep));
    SolverTaskDescriptionNode.appendChild(TimeStepNode);
    
    ErrorToleranceNode = docNode.createElement('ErrorTolerance');
    ErrorToleranceNode.setAttribute('Absolut',num2str(VCMLAbsoluteError));
    ErrorToleranceNode.setAttribute('Relative',num2str(VCMLRelativeError));
    SolverTaskDescriptionNode.appendChild(ErrorToleranceNode);
    
    OutputOptionsNode = docNode.createElement('OutputOptions');
    OutputOptionsNode.setAttribute('OutputTimeStep',num2str(VCMLOutputTimeStep));
    SolverTaskDescriptionNode.appendChild(OutputOptionsNode);
    
    SundialsSolverOptionsNode = docNode.createElement('SundialsSolverOptions');
    SolverTaskDescriptionNode.appendChild(SundialsSolverOptionsNode);
    
    maxOrderAdvectionNode = docNode.createElement('maxOrderAdvection');
    SundialsSolverOptionsNode.appendChild(maxOrderAdvectionNode);
    maxOrderAdvectionNode.setTextContent('2');
    
    NumberProcessorsNode = docNode.createElement('NumberProcessors');
    SolverTaskDescriptionNode.appendChild(NumberProcessorsNode);
    NumberProcessorsNode.setTextContent(num2str(1));
    
    MathOverridesNode = docNode.createElement('MathOverrides');
    SimulationNode.appendChild(MathOverridesNode);
    
    MeshSpecificationNode = docNode.createElement('MeshSpecification');
    SimulationNode.appendChild(MeshSpecificationNode);
    
    SizeNode = docNode.createElement('Size');
    SizeNode.setAttribute('X',num2str(size(all_compartments_image,2)));
    SizeNode.setAttribute('Y',num2str(size(all_compartments_image,1)));
    SizeNode.setAttribute('Z',num2str(size(all_compartments_image,3)));
    MeshSpecificationNode.appendChild(SizeNode);
    
    
    % SpatialObjects
    
    removeChildrenByTagName(SimulationSpecNode, 'SpatialObjects');
    
    SpatialObjectsNode = docNode.createElement('SpatialObjects');
    SimulationSpecNode.appendChild(SpatialObjectsNode);
    
    for j = 1:length(all_compartments_objects)
        name = all_compartments_objects_names{j};
        object = all_compartments_objects(name);
        
        SpatialObjectNode = docNode.createElement('SpatialObject');
        SpatialObjectNode.setAttribute('Name',object.SpatialObjectName);
        SpatialObjectNode.setAttribute('Type','Volume');
        SpatialObjectNode.setAttribute('subVolume',object.subVolumeName);
        SpatialObjectNode.setAttribute('regionId',num2str(j-1));
        SpatialObjectsNode.appendChild(SpatialObjectNode);
        
        QuantityCategoryListNode = docNode.createElement('QuantityCategoryList');
        SpatialObjectNode.appendChild(QuantityCategoryListNode);
        
        QuantityCategoryNode = docNode.createElement('QuantityCategory');
        QuantityCategoryNode.setAttribute('Name','VolumeCentroid');
        QuantityCategoryNode.setAttribute('Enabled',mat2str(false));
        QuantityCategoryListNode.appendChild(QuantityCategoryNode);
        
        QuantityCategoryNode = docNode.createElement('QuantityCategory');
        QuantityCategoryNode.setAttribute('Name','InteriorVelocity');
        QuantityCategoryNode.setAttribute('Enabled',mat2str(false));
        QuantityCategoryListNode.appendChild(QuantityCategoryNode);
        
        QuantityCategoryNode = docNode.createElement('QuantityCategory');
        QuantityCategoryNode.setAttribute('Name','VolumeRegionSize');
        QuantityCategoryNode.setAttribute('Enabled',mat2str(true));
        QuantityCategoryListNode.appendChild(QuantityCategoryNode);
    end
    
    for j = 1:length(all_object_membrane_VCML_data)
        VCML_data = all_object_membrane_VCML_data(j);
        
        SpatialObjectNode = docNode.createElement('SpatialObject');
        SpatialObjectNode.setAttribute('Name',VCML_data.SpatialObjectName);
        SpatialObjectNode.setAttribute('Type','Surface');
        SpatialObjectNode.setAttribute('subVolumeOutside',VCML_data.subVolumeOutside);
        SpatialObjectNode.setAttribute('subVolumeInside',VCML_data.subVolumeInside);
        SpatialObjectNode.setAttribute('regionIdOutside',num2str(VCML_data.regionIdOutside));
        SpatialObjectNode.setAttribute('regionIdInside',num2str(VCML_data.regionIdInside));
        SpatialObjectsNode.appendChild(SpatialObjectNode);
        
        QuantityCategoryListNode = docNode.createElement('QuantityCategoryList');
        SpatialObjectNode.appendChild(QuantityCategoryListNode);
        
        QuantityCategoryNode = docNode.createElement('QuantityCategory');
        QuantityCategoryNode.setAttribute('Name','SurfaceNormal');
        QuantityCategoryNode.setAttribute('Enabled',mat2str(false));
        QuantityCategoryListNode.appendChild(QuantityCategoryNode);
        
        QuantityCategoryNode = docNode.createElement('QuantityCategory');
        QuantityCategoryNode.setAttribute('Name','SurfaceVelocity');
        QuantityCategoryNode.setAttribute('Enabled',mat2str(false));
        QuantityCategoryListNode.appendChild(QuantityCategoryNode);
        
        QuantityCategoryNode = docNode.createElement('QuantityCategory');
        QuantityCategoryNode.setAttribute('Name','DistanceToSurface');
        QuantityCategoryNode.setAttribute('Enabled',mat2str(false));
        QuantityCategoryListNode.appendChild(QuantityCategoryNode);
        
        QuantityCategoryNode = docNode.createElement('QuantityCategory');
        QuantityCategoryNode.setAttribute('Name','DirectionToSurface');
        QuantityCategoryNode.setAttribute('Enabled',mat2str(false));
        QuantityCategoryListNode.appendChild(QuantityCategoryNode);
        
        QuantityCategoryNode = docNode.createElement('QuantityCategory');
        QuantityCategoryNode.setAttribute('Name','SurfaceSize');
        QuantityCategoryNode.setAttribute('Enabled',mat2str(true));
        QuantityCategoryListNode.appendChild(QuantityCategoryNode);
    end
    
    
    % MicroscopeMeasurement
    
    removeChildrenByTagName(SimulationSpecNode, 'MicroscopeMeasurement');
    
    MicroscopeMeasurementNode = docNode.createElement('MicroscopeMeasurement');
    MicroscopeMeasurementNode.setAttribute('Name','fluor');
    SimulationSpecNode.appendChild(MicroscopeMeasurementNode);
    
    ConvolutionKernelNode = docNode.createElement('ConvolutionKernel');
    ConvolutionKernelNode.setAttribute('Type','ProjectionZKernel');
    MicroscopeMeasurementNode.appendChild(ConvolutionKernelNode);
    
    % ConvolutionKernelNode.setAttribute('Type','GaussianConvolutionKernel');
    % SigmaXYNode = docNode.createElement('SigmaXY');
    % ConvolutionKernelNode.appendChild(SigmaXYNode);
    % SigmaXYNode.setTextContent(num2str(0));
    % SigmaZNode = docNode.createElement('SigmaZ');
    % ConvolutionKernelNode.appendChild(SigmaZNode);
    % SigmaZNode.setTextContent(num2str(0));
    % MicroscopeMeasurementNode.appendChild(ConvolutionKernelNode);
    
    % relationshipModelNode = docNode.createElement('relationshipModel');
    % BioModelNode.appendChild(relationshipModelNode);
    
    % vcmetadataNode = docNode.createElement('vcmetadata');
    % BioModelNode.appendChild(vcmetadataNode);
end


warning('CellOrganizer:instance2VCML', 'TODO: Double check unit adjustments (um_to_vcml_length_unit, etc.)');


%%%
%%%Save the file

xmlwrite(savepath, docNode);
zip([savepath, '.zip'], savepath);

result = true;

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
% Find adjacent values using NET file info

adjacent_pairs = zeros(0, 2);
network_info_compartments_names = network_info_compartments.keys();
for i = 1:length(network_info_compartments)
    network_info_compartment_name = network_info_compartments_names{i};
    network_info_compartment = network_info_compartments(network_info_compartment_name);
    network_info_compartment_spatial_dimensions = network_info_compartment.spatial_dimensions;
    
    if network_info_compartment_spatial_dimensions ~= 3
        continue;
    end
    
    if ~isempty(network_info_compartment.outside) && network_info_compartments.isKey(network_info_compartment.outside)
        network_info_compartment_actual_outside = network_info_compartment.outside;
        
        network_info_compartment_actual_outside = network_info_compartments(network_info_compartment_actual_outside).outside;
        
        if ~isempty(network_info_compartment_actual_outside)
            if all_compartment_indices.isKey(network_info_compartment_name) && all_compartment_indices.isKey(network_info_compartment_actual_outside)
                    adjacent_pairs(end+1, :) = [all_compartment_indices(network_info_compartment_name), all_compartment_indices(network_info_compartment_actual_outside)];
            end
        end
    end
end
adjacent_pairs = sortrows(adjacent_pairs);
adjacent_pairs = sort(adjacent_pairs, 2);
adjacent_pairs = unique(adjacent_pairs, 'rows');
adjacent_pairs = adjacent_pairs(adjacent_pairs(:, 1) ~= adjacent_pairs(:, 2), :);
end




function value = translateWithDefaultIdentity(map_object, key)

if map_object.isKey(key)
    value = map_object(key);
else
    value = key;
end

end




function value = inCellArray(cell_array, key)

value = any(cellfun(@(x)valuesEqual(x, key), cell_array));

end




function value = intersectCellArrays(cell_array1, cell_array2)

value = {};
for i = 1:length(cell_array1)
    if inCellArray(cell_array1{i}, cell_array2)
        value{end+1} = cell_array1{i};
    end
end

end




function removeChildrenByTagName(node, tag_name)

children = node.getElementsByTagName(tag_name);
for i = children.getLength():-1:1
    child = children.item(i-1);
    child.getParentNode().removeChild(child);
end

end




function result = checkForChildrenByTagName(node, tag_name)

result = node.getElementsByTagName(tag_name).getLength() > 0;

end




function result = getFirstChildByTagName(node, tag_name)

result = node.getElementsByTagName(tag_name).item(0);

end




function result = getChild(node, tag_name, match_attributes_struct)

if nargin < 3 || isempty(match_attributes_struct)
    match_attributes_struct = struct();
end

result = [];

match_attributes_names = fieldnames(match_attributes_struct);
elements = node.getElementsByTagName(tag_name);
for i = 1:elements.getLength()
    element = elements.item(i-1);
    is_match = true;
    for j = 1:length(match_attributes_names)
        match_attribute_name = match_attributes_names{j};
        match_attribute = match_attributes_struct.(match_attribute_name);
        element_match_attribute = getAttributeValue(element, match_attribute_name, class(match_attribute));
        if ~valuesEqual(match_attribute, element_match_attribute)
            is_match = false;
            break;
        end
    end
    if is_match
        result = element;
        break;
    end
end


end




function result = createOrGetChild(node, tag_name, match_attributes_struct, set_attributes_struct)

if nargin < 3 || isempty(match_attributes_struct)
    match_attributes_struct = struct();
end
if nargin < 4 || isempty(set_attributes_struct)
    set_attributes_struct = struct();
end

result = getChild(node, tag_name, match_attributes_struct);

if isempty(result)
    result = node.getOwnerDocument().createElement(tag_name);
    node.appendChild(result);

    match_attributes_names = fieldnames(match_attributes_struct);
    for i = 1:length(match_attributes_names)
        match_attribute_name = match_attributes_names{i};
        match_attribute = match_attributes_struct.(match_attribute_name);
        result.setAttribute(match_attribute_name, match_attribute);
    end
end

set_attributes_names = fieldnames(set_attributes_struct);
for i = 1:length(set_attributes_names)
    set_attribute_name = set_attributes_names{i};
    set_attribute = set_attributes_struct.(set_attribute_name);
    result.setAttribute(set_attribute_name, set_attribute);
end


end




function attribute = javaStringToValue(attribute, type)

attribute = char(attribute);
if any(strcmpi(type, {'str', 'string', 'char'}))
    % Do nothing
elseif any(strcmpi(type, {'bool', 'boolean'}))
    attribute = strcmpi(attribute, 'true');
elseif any(strcmpi(type, {'num', 'number', 'int', 'integer', 'float', 'single', 'double'}))
    attribute = str2double(attribute);
else
    error('Unknown type ''%s''', type);
end

end




function attribute = getAttributeValue(node, attribute_name, type)

attribute = javaStringToValue(node.getAttribute(attribute_name), type);

end




function attribute = getTextContentValue(node, type)

attribute = javaStringToValue(node.getTextContent(), type);

end




function result = valuesEqual(value1, value2)

if any(strcmpi(class(value1), {'str', 'string', 'char'}))
    result = strcmp(value1, value2);
else
    result = value1 == value2;
end

end




function network_info = readNetwork(NETfile)

parameters = containers.Map('KeyType', 'char', 'ValueType', 'any');
compartments = containers.Map('KeyType', 'char', 'ValueType', 'any');
molecule_types = struct('index', [], 'type', {}, 'comment', {});
observables = struct('index', [], 'type', {}, 'name', {}, 'observable_expressions', {}, 'comment', {});
species = struct('index', [], 'species_graph', {}, 'concentration', {}, 'compartments', {}, 'diffusion_coefficient', [], 'comment', {});
reaction_rules = struct('name', {}, 'reaction', {}, 'direction', {}, 'forward_rate_constant', {}, 'reverse_rate_constant', {}, 'comment', {});
reactions = struct('index', {}, 'reactant_indices', {}, 'product_indices', {}, 'rate_constant', {}, 'rate_constant_parameters', {}, 'unit_conversion', {}, 'comment', {});

if isempty(NETfile) || ~ischar(NETfile)
    % Produce skeleton anyway
    contents = '';
else
    contents = fileread(NETfile);
end
contents_lines = strsplit(contents, {'\n', '\r'});
section = '';

warning('CellOrganizer:instance2VCML', 'Check ''reverse_rate_constant''')


bng_patterns = getBNGRegexps();

section_patterns_options = struct('begins_line', true, 'ends_line', true);
reaction_rules_pattern_options = struct('begins_line', true, 'ends_line', true, 'separator', '  ', 'separator_has_variable_length', false);

section_patterns = containers.Map('KeyType', 'char', 'ValueType', 'char');
section_patterns('parameters') = createRegexp(struct('name', {'index', 'name', 'value_expression'}, 'type', {'integer_signless', 'identifier', 'expression_spaceless'}, 'required', num2cell(logical([0, 1, 1]))), section_patterns_options);
section_patterns('compartments') = createRegexp(struct('name', {'name', 'spatial_dimensions', 'size_expression', 'outside'}, 'type', {'identifier', 'integer', 'expression_spaceless', 'identifier'}, 'required', num2cell(logical([1, 1, 0, 0]))), section_patterns_options);
section_patterns('molecule types') = createRegexp(struct('name', {'index', 'type'}, 'type', {'integer_signless', 'custom'}, 'pattern', {'', bng_patterns.net_molecule_type}, 'required', num2cell(logical([0, 1]))), section_patterns_options);
section_patterns('observables') = createRegexp(struct('name', {'index', 'type', 'name', 'observable_expressions'}, 'type', {'integer_signless', 'identifier', 'identifier', 'custom'}, 'pattern', {'', '', '', bng_patterns.net_observables_pattern_list}, 'required', num2cell(logical([0, 1, 1, 1]))), section_patterns_options);
section_patterns('species') = createRegexp(struct('name', {'index', 'species_graph', 'concentration'}, 'type', {'integer_signless', 'custom', 'expression_spaceless'}, 'pattern', {'', bng_patterns.net_pre_species_def, ''}, 'required', num2cell(logical([0, 1, 1]))), section_patterns_options);
section_patterns('reaction rules') = createRegexp(struct('name', {'name', 'reaction', 'rate_constants'}, 'type', {'custom', 'custom', 'custom'}, 'pattern', {bng_patterns.net_reaction_rules_name, bng_patterns.net_reaction_rules_rule, bng_patterns.net_reaction_rules_rate_constant}, 'required', num2cell(logical([0, 1, 1]))), reaction_rules_pattern_options);
section_patterns('reactions') = createRegexp(struct('name', {'index', 'reactant_indices', 'product_indices', 'rate_constant'}, 'type', {'custom', 'custom', 'custom', 'custom'}, 'pattern', {bng_patterns.net_reaction_index, bng_patterns.net_reaction_reactant_indices, bng_patterns.net_reaction_product_indices, bng_patterns.net_reaction_rate_constant}, 'required', num2cell(logical([1, 1, 1, 1]))), section_patterns_options);

% section_patterns('groups') = createRegexp(struct());

species_compartments_pattern = ['@' bng_patterns.string_pattern ''];
rate_constant_parameters_pattern = bng_patterns.string_pattern;
value_expression_parameters_pattern = rate_constant_parameters_pattern;
species_comment_diffusion_coefficient_pattern = ['diffusion_coefficient=(' bng_patterns.float_pattern ')'];
reaction_comment_unit_conversion_pattern = ['unit_conversion=(' bng_patterns.net_expression ')'];

default_diffusion_coefficient = nan;

for i = 1:length(contents_lines)
    contents_line = strtrim(contents_lines{i});
    contents_line_and_comment = strsplit(contents_line, '#');
    if length(contents_line_and_comment) > 1
        contents_line = contents_line_and_comment{1};
        contents_line_comment = strjoin(contents_line_and_comment(2:end), '#');
        contents_line_comment = strtrim(contents_line_comment);
        contents_line = strtrim(contents_line);
    else
        contents_line_comment = '';
    end
    clear contents_line_and_comment;
    
    % Comments, blank lines, and section delimiters
    if length(contents_line) == 0
        continue;
    elseif strcmp(contents_line, 'begin parameters')
        section = 'parameters';
        continue;
    elseif strcmp(contents_line, 'end parameters')
        section = '';
        continue;
    elseif strcmp(contents_line, 'begin compartments')
        section = 'compartments';
        continue;
    elseif strcmp(contents_line, 'end compartments')
        section = '';
        continue;
    elseif strcmp(contents_line, 'begin molecule types')
        section = 'molecule types';
        continue;
    elseif strcmp(contents_line, 'end molecule types')
        section = '';
        continue;
    elseif strcmp(contents_line, 'begin observables')
        section = 'observables';
        continue;
    elseif strcmp(contents_line, 'end observables')
        section = '';
        continue;
    elseif strcmp(contents_line, 'begin energy patterns')
        section = 'energy patterns';
        warning('CellOrganizer:instance2VCML', 'Energy patterns ignored');
        continue;
    elseif strcmp(contents_line, 'end energy patterns')
        section = '';
        continue;
    elseif strcmp(contents_line, 'begin functions')
        section = 'functions';
        warning('CellOrganizer:instance2VCML', 'Functions ignored');
        continue;
    elseif strcmp(contents_line, 'end functions')
        section = '';
        continue;
    elseif any(strcmp(contents_line, {'begin species', 'begin seed species'}))
        section = 'species';
        % Diffusion coefficients have to be added manually for now
        default_diffusion_coefficient_tokens = regexp(contents_line_comment, species_comment_diffusion_coefficient_pattern, 'tokens');
        if length(default_diffusion_coefficient_tokens) > 0
            default_diffusion_coefficient = str2double(default_diffusion_coefficient_tokens{1}{1});
        end
        continue;
    elseif any(strcmp(contents_line, {'end species', 'end seed species'}))
        section = '';
        continue;
    elseif strcmp(contents_line, 'begin reaction rules')
        section = 'reaction rules';
        continue;
    elseif strcmp(contents_line, 'end reaction rules')
        section = '';
        continue;
    elseif strcmp(contents_line, 'begin reactions')
        section = 'reactions';
        continue;
    elseif strcmp(contents_line, 'end reactions')
        section = '';
        continue;
    elseif strcmp(contents_line, 'begin groups')
        section = 'groups';
        continue;
    elseif strcmp(contents_line, 'end groups')
        section = '';
        continue;
    end
    
    switch section
        case 'parameters'
            info_struct = regexp(contents_line, section_patterns('parameters'), 'names');
            info_struct.index = str2double(info_struct.index);
            info_struct.name = strtrim(info_struct.name);
            info_struct.value_expression_parameters = regexp(info_struct.value_expression, value_expression_parameters_pattern, 'match');
            
            info_struct.comment = contents_line_comment;
            parameters(info_struct.name) = info_struct;
            
        case 'compartments'
            info_struct = regexp(contents_line, section_patterns('compartments'), 'names');
            info_struct.name = strtrim(info_struct.name);
            info_struct.spatial_dimensions = str2double(info_struct.spatial_dimensions);
            info_struct.size_expression = strtrim(info_struct.size_expression);
            info_struct.outside = strtrim(info_struct.outside);
            
            info_struct.comment = contents_line_comment;
            compartments(info_struct.name) = info_struct;
            
        case 'molecule types'
            info_struct = regexp(contents_line, section_patterns('molecule types'), 'names');
            info_struct.index = str2double(info_struct.index);
            if isnan(info_struct.index)
                info_struct.index = length(molecule_types)+1;
            end
            info_struct.type = strtrim(info_struct.type);
            
            info_struct.comment = contents_line_comment;
            molecule_types(end+1) = info_struct;
            
        case 'observables'
            % Do nothing
            % warning('CellOrganizer:instance2VCML', 'Unfinished');
            %{
            info_struct = regexp(contents_line, section_patterns('observables'), 'names');
            info_struct.index = str2double(info_struct.index);
            if isnan(info_struct.index)
                info_struct.index = length(observables)+1;
            end
            info_struct.name = strtrim(info_struct.name);
            info_struct.observable_expressions = strsplit(strtrim(info_struct.observable_expressions));
            
            info_struct.comment = contents_line_comment;
            observables(end+1) = info_struct;
            %}
            
        case 'species'
            info_struct = regexp(contents_line, section_patterns('species'), 'names');
            info_struct.index = str2double(info_struct.index);
            if isnan(info_struct.index)
                info_struct.index = length(species)+1;
            end
            info_struct.species_graph = strtrim(info_struct.species_graph);
            info_struct.concentration = strtrim(info_struct.concentration);
            info_struct.compartments = regexp(info_struct.species_graph, species_compartments_pattern, 'match');
            info_struct.compartments = strrep(info_struct.compartments, '@', '');
            info_struct.compartments = unique(info_struct.compartments);
            % Diffusion coefficients have to be added manually for now
            diffusion_coefficient_tokens = regexp(contents_line_comment, species_comment_diffusion_coefficient_pattern, 'tokens');
            if length(diffusion_coefficient_tokens) > 0
                info_struct.diffusion_coefficient = str2double(diffusion_coefficient_tokens{1}{1});
            else
                info_struct.diffusion_coefficient = default_diffusion_coefficient;
            end
            
            info_struct.comment = contents_line_comment;
            species(end+1) = info_struct;
            
        case 'reaction rules'
            % Do nothing
            % warning('CellOrganizer:instance2VCML', 'Unfinished');
            %{
            info_struct = regexp(contents_line, section_patterns('reaction rules'), 'names');
            info_struct.name = strtrim(info_struct.name);
            info_struct.reaction = strtrim(info_struct.reaction);
            if length(strfind(info_struct.reaction, '<->')) > 0
                info_struct.direction = 'both';
            else
                info_struct.direction = 'forward';
            end
            info_struct.rate_constants = strtrim(strsplit(info_struct.rate_constants, ','));
            info_struct.forward_rate_constant = info_struct.rate_constants{1};
            if length(info_struct.rate_constants) > 1
                info_struct.reverse_rate_constant = info_struct.rate_constants{2};
            elseif strcmp(info_struct.direction, 'forward')
                info_struct.reverse_rate_constant = '0';
            else
                % TODO: Correct this. How does BioNetGen compute the reverse rate_constant?
                info_struct.reverse_rate_constant = info_struct.forward_rate_constant;
            end
            info_struct = rmfield(info_struct, 'rate_constants');
            
            info_struct.comment = contents_line_comment;
            reaction_rules(end+1) = info_struct;
            %}
            
        case 'reactions'
            info_struct = regexp(contents_line, section_patterns('reactions'), 'names');
            info_struct.index = str2double(info_struct.index);
            if isnan(info_struct.index)
                info_struct.index = length(reactions)+1;
            end
            info_struct.reactant_indices = str2double(strsplit(info_struct.reactant_indices, ','));
            info_struct.product_indices = str2double(strsplit(info_struct.product_indices, ','));
            % 0 means nothing on that side
            info_struct.reactant_indices = info_struct.reactant_indices(info_struct.reactant_indices > 0);
            info_struct.product_indices = info_struct.product_indices(info_struct.product_indices > 0);
            info_struct.rate_constant_parameters = regexp(info_struct.rate_constant, rate_constant_parameters_pattern, 'match');
            info_struct.unit_conversion = regexp(contents_line_comment, reaction_comment_unit_conversion_pattern, 'tokens');
            if length(info_struct.unit_conversion) > 0
                info_struct.unit_conversion = info_struct.unit_conversion{1}{1};
            else
                info_struct.unit_conversion = '1';
            end
            
            info_struct.comment = contents_line_comment;
            reactions(end+1) = info_struct;
            
        case 'groups'
            % Do nothing
            % warning('CellOrganizer:instance2VCML', 'Unfinished');
            %{
            info_struct = regexp(contents_line, section_patterns('groups'), 'names');
            % info_struct.index = str2double(info_struct.index);
            % info_struct.observable_expressions = strsplit(observable_expressions);
            
            info_struct.comment = contents_line_comment;
            observables(info_struct.name) = info_struct;
            %}
            
        case 'energy patterns'
            % Do nothing
        case 'functions'
            % Do nothing
            
        otherwise
            error('Network file unreadable');
    end
end

if length(contents) > 0 && length(compartments) == 0
    [NETfile_dir, NETfile_name, NETfile_ext] = fileparts(NETfile);
    error(sprintf(['\n', ...
        '\nNo or empty compartments section in NET file.', ...
        '\nPlease generate NET files through BNG2.pl with commands similar to the following:', ...
        '\n', ...
        '\n    load %s.bngl', ...
        '\n    action generate_network({overwrite=>1})', ...
        '\n    action writeFile({format=>\"net\", overwrite=>1, pretty_formatting=>1})', ...
        '\n'], NETfile_name));
end

network_info = struct();
network_info.parameters = parameters;
network_info.compartments = compartments;
network_info.molecule_types = molecule_types;
% network_info.observables = observables;
network_info.species = species;
% network_info.reaction_rules = reaction_rules;
network_info.reactions = reactions;
% network_info.groups = groups;

network_info.regexp_patterns = struct();
network_info.regexp_patterns.bng_patterns = bng_patterns;
network_info.regexp_patterns.section_patterns = section_patterns;
network_info.regexp_patterns.species_compartments_pattern = species_compartments_pattern;
network_info.regexp_patterns.rate_constant_parameters_pattern = rate_constant_parameters_pattern;
network_info.regexp_patterns.value_expression_parameters_pattern = value_expression_parameters_pattern;
network_info.regexp_patterns.reaction_comment_unit_conversion_pattern = reaction_comment_unit_conversion_pattern;

end
