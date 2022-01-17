function [ result ] = instance2MCellMDL( CSGdata, meshData, models, imgs, options, savepath )
%INSTANCE2MCELLMDL Writes geometry and compartmental reactions to MCellMDL file.
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

% 2019-06-11 Taraz Buck - Copied `instance2VCML.m` to `instance2MCellMDL.m`.


warning('off', 'CellOrganizer:instance2MCellMDL');

result = false;

if nargin < 3
    % warning('CellOrganizer:instance2MCellMDL:missingRequiredArgument', 'First three arguments are required');
    warning('CellOrganizer:instance2MCellMDL', 'First three arguments are required');
    return;
end

if nargin < 4
    % warning('CellOrganizer:instance2MCellMDL:missingOptionalArgument', 'Argument savepath not given, defaulting to ''./model.mdl''');
    warning('CellOrganizer:instance2MCellMDL', 'Argument savepath not given, defaulting to ''./model.mdl''');
    savepath = './model.mdl';
end
[savepath_dir, savepath_name, savepath_ext] = fileparts(savepath);
savepath_dir = strrep(savepath_dir, pwd(), '.');
savepath = fullfile(savepath_dir, [savepath_name, savepath_ext]);
script_savepath = fullfile(savepath_dir, [savepath_name, '.sh']);

imgs_valid = true(length(imgs), 1);
for i = 1:length(imgs)
    imgs_valid(i) = imgs_valid(i) && numel(imgs{i}) > 0;
    imgs_valid(i) = imgs_valid(i) && ndims(imgs{i}) == 3;
end
if ~all(imgs_valid)
    % warning('CellOrganizer:instance2MCellMDL:imgsEmpty', 'Argument imgs must contain nonempty 3D arrays');
    warning('CellOrganizer:instance2MCellMDL', 'Argument imgs must contain nonempty 3D arrays');
    return;
end


mcell_length_unit = UnitsManager('um');
mcell_area_unit = mcell_length_unit^2;
mcell_volume_unit = mcell_length_unit^3;
mcell_concentration_unit = UnitsManager('M');
mcell_density_unit = UnitsManager('molecule.um-2');
mcell_time_unit = UnitsManager('s');
mcell_microscopic_reaction_rate_unit = UnitsManager('molecule.s-1');
mcell_volume_reaction_rate_unit = mcell_microscopic_reaction_rate_unit;
mcell_reaction_rate_base_unit = mcell_microscopic_reaction_rate_unit;
mcell_membrane_reaction_rate_variable_unit = mcell_density_unit^-1;
mcell_volume_reaction_rate_variable_unit = mcell_concentration_unit^-1;
mcell_reaction_rate_unit_volume_exponent_function = @(given_n_reactants, given_n_reactants_at_surfaces)max(0, given_n_reactants - max(0, given_n_reactants_at_surfaces - 1) - 1);
mcell_reaction_rate_unit_volume_part_function = @(given_n_reactants, given_n_reactants_at_surfaces)mcell_volume_reaction_rate_variable_unit^mcell_reaction_rate_unit_volume_exponent_function(given_n_reactants, given_n_reactants_at_surfaces);
mcell_reaction_rate_unit_membrane_exponent_function = @(given_n_reactants, given_n_reactants_at_surfaces)max(0, given_n_reactants_at_surfaces - 1);
mcell_reaction_rate_unit_membrane_part_function = @(given_n_reactants, given_n_reactants_at_surfaces)mcell_membrane_reaction_rate_variable_unit^mcell_reaction_rate_unit_membrane_exponent_function(given_n_reactants, given_n_reactants_at_surfaces);
mcell_reaction_rate_unit_function = @(given_n_reactants, given_n_reactants_at_surfaces)mcell_reaction_rate_base_unit * mcell_reaction_rate_unit_volume_part_function(given_n_reactants, given_n_reactants_at_surfaces) * mcell_reaction_rate_unit_membrane_part_function(given_n_reactants, given_n_reactants_at_surfaces) * given_n_reactants^(-max(0, given_n_reactants_at_surfaces - 1));
mcell_reaction_rate_constant_unit_function = @(given_n_reactants, given_n_reactants_at_surfaces)UnitsManager('M')^(-mcell_reaction_rate_unit_volume_exponent_function(given_n_reactants, given_n_reactants_at_surfaces)) * (mcell_density_unit * given_n_reactants)^-mcell_reaction_rate_unit_membrane_exponent_function(given_n_reactants, given_n_reactants_at_surfaces) / UnitsManager('s');
mcell_diffusion_coefficient_unit = UnitsManager('cm2.s-1');

n_partitions = 20;
index_to_axis = {'X', 'Y', 'Z'};
axis_to_index = containers.Map({'X', 'Y', 'Z'}, [1, 2, 3]);



% Process NET file and combine with geometry

options.output.NET.output_length_unit = mcell_length_unit;
options.output.NET.output_concentration_unit = mcell_concentration_unit;
options.output.NET.output_time_unit = mcell_time_unit;

options.output.NET.use_image_adjacency = false;
options.output.temp = struct();
options.output.temp.savepath = savepath;

options = nameModels(CSGdata, meshData, models, imgs, options);
options.output.use_individual_meshes = true;
[CSGdata, meshData] = convertCSGToMesh(CSGdata, meshData, models, imgs, options);

network_with_geometry_info = readNetworkIntoGeometry(CSGdata, meshData, models, imgs, options);
network_info = network_with_geometry_info.network_info;
geometry_info = network_with_geometry_info.geometry_info;
meshData = geometry_info.meshData;

% network_with_geometry_info
% network_info
% geometry_info
% meshData

models = geometry_info.models;
imgs = geometry_info.imgs;
meshes = geometry_info.meshes;
dim_chars = geometry_info.dim_chars;
sign_chars = geometry_info.sign_chars;
connectivity = geometry_info.connectivity;
connectivity_se = geometry_info.connectivity_se;
avogadro_constant = geometry_info.avogadro_constant;
avogadro_constant_value_expression = geometry_info.avogadro_constant_value_expression;
avogadro_constant_units = geometry_info.avogadro_constant_units;
model_names = geometry_info.model_names;
model_cytonuclearflags = geometry_info.model_cytonuclearflags;
name_map = geometry_info.name_map;
named_imgs = geometry_info.named_imgs;
named_imgs_cytonuclearflags = geometry_info.named_imgs_cytonuclearflags;
names = geometry_info.names;
all_compartments_image_before_downsampling = geometry_info.all_compartments_image_before_downsampling;
all_compartments_image = geometry_info.all_compartments_image;
named_imgs_before_downsampling = geometry_info.named_imgs_before_downsampling;
adjacent_pairs = geometry_info.adjacent_pairs;
n_compartments = geometry_info.n_compartments;
all_compartment_indices = geometry_info.all_compartment_indices;
n_network_info_compartments = geometry_info.n_network_info_compartments;
network_info_compartments_keys = geometry_info.network_info_compartments_keys;
network_info_compartments_keys_indices_map = geometry_info.network_info_compartments_keys_indices_map;
network_info_adjacency_max_degree = geometry_info.network_info_adjacency_max_degree;
network_info_adjacency_matrix = geometry_info.network_info_adjacency_matrix;
network_info_adjacency_matrices = geometry_info.network_info_adjacency_matrices;
is_network_info_compartment_pair_adjacent = geometry_info.is_network_info_compartment_pair_adjacent;
adjacent_pairs_max_degree = geometry_info.adjacent_pairs_max_degree;
adjacent_pairs_matrix = geometry_info.adjacent_pairs_matrix;
is_compartment_pair_adjacent = geometry_info.is_compartment_pair_adjacent;
all_compartments_objects = geometry_info.all_compartments_objects;
all_compartments_objects_image = geometry_info.all_compartments_objects_image;
all_compartments_objects_names = geometry_info.all_compartments_objects_names;
all_compartments_objects_names_to_compartment_names = geometry_info.all_compartments_objects_names_to_compartment_names;
all_compartments_objects_indices = geometry_info.all_compartments_objects_indices;
n_compartments_objects = geometry_info.n_compartments_objects;
object_adjacent_pairs = geometry_info.object_adjacent_pairs;
getAdjacentValues = geometry_info.getAdjacentValues;
getAdjacentCompartments = geometry_info.getAdjacentCompartments;

all_compartment_data = geometry_info.all_compartment_data;
all_object_data = geometry_info.all_object_data;
all_object_membrane_data = geometry_info.all_object_membrane_data;
all_membrane_data = geometry_info.all_membrane_data;
all_membrane_names = geometry_info.all_membrane_names;
is_membrane_function = geometry_info.is_membrane_function;
species_index_to_name_function = geometry_info.species_index_to_name_function;

translateWithDefaultIdentity = geometry_info.translateWithDefaultIdentity;
parameterReferenceBase = geometry_info.parameterReferenceBase;
evaluateExpression = geometry_info.evaluateExpression;


function result = parameterReference(given_parameter_name)
    result = parameterReferenceBase(given_parameter_name, network_info.parameters);
end



NETfile = options.output.NET.filename;

NETTranslations = options.output.NET.translations;
NETDownsample = options.output.NET.downsampling;
NETAddTranslocationIntermediates = options.output.NET.add_translocation_intermediates;
NETDefaultDiffusionCoefficient = options.output.NET.default_diffusion_coefficient;

MCellMDLNumSimulations = options.output.MCellMDL.num_simulations;
MCellMDLEndTime = options.output.MCellMDL.end_time;
MCellMDLDefaultTimeStep = options.output.MCellMDL.default_time_step;
MCellMDLMaxTimeStep = options.output.MCellMDL.max_time_step;
MCellMDLOutputTimeStep = options.output.MCellMDL.output_time_step;
MCellMDLInteractionRadius = options.output.MCellMDL.interaction_radius;
MCellMDLMaxRealTime = options.output.MCellMDL.max_real_time;
MCellMDLUseReactionRateHack = options.output.MCellMDL.use_reaction_rate_hack;
MCellMDLInputFilenamePattern = options.output.MCellMDL.input_filename_pattern;

net_length_unit = options.output.NET.units.length;
net_area_unit = unitsPower(net_length_unit, 2);
net_volume_unit = unitsPower(net_length_unit, 3);
net_time_unit = options.output.NET.units.time;
net_concentration_unit = options.output.NET.units.concentration;
net_volume_substance_unit = [net_concentration_unit,'.',net_volume_unit];
net_membrane_substance_unit = 'molecule';
net_diffusion_coefficient_unit = [net_area_unit, '.', unitsPower(net_time_unit, -1)];

% warning('CellOrganizer:instance2MCellMDL:assumingUnits', 'Assuming NET file concentrations are extensive');
warning('CellOrganizer:instance2MCellMDL', 'Assuming NET file concentrations are extensive');

% avogadro_constant = 6.022140857e23;

function result2 = mcell_reaction_rate_info_function(given_reaction)
    given_n_reactants = length(given_reaction.reactant_indices);
    given_n_reactants_at_surfaces = length(given_reaction.reactants_indices_at_surfaces);
    given_n_reactants_in_volumes = given_n_reactants - given_n_reactants_at_surfaces;
    given_compartment = given_reaction.compartment;
    given_compartment_volume = network_info.compartments(given_compartment).size_expression;
    result2 = struct();
    result2.concentration_exponent = 0;
    result2.area_exponent = 0;
    result2.N_Avo_exponent = 0;
    % result2.N_Avo_exponent = given_n_reactants_at_surfaces;
    result2.n_reactants_exponent = 0;
    result2.volume_exponent = mcell_reaction_rate_unit_volume_exponent_function(given_n_reactants, given_n_reactants_at_surfaces);
    result2.membrane_exponent = mcell_reaction_rate_unit_membrane_exponent_function(given_n_reactants, given_n_reactants_at_surfaces);
    result2.concentration_exponent = -result2.volume_exponent;
    result2.area_exponent = result2.membrane_exponent;
    result2.n_reactants_exponent = -result2.membrane_exponent;
    result2.N_Avo_exponent = -result2.volume_exponent;
    
    result2.rate_scale = DimensionedExpression(1);
    result2.rate_scale = result2.rate_scale * parameterReference('N_Avo')^result2.N_Avo_exponent * DimensionedExpression(given_n_reactants)^result2.n_reactants_exponent;
    result2.rate_scale = result2.rate_scale * network_info.compartments(given_reaction.compartment).size_expression^(given_n_reactants - 1);
    result2.rate_scale = result2.rate_scale / parameterReference('eff_width')^result2.area_exponent;
    
    result2.rate_units = mcell_reaction_rate_unit_function(given_n_reactants, given_n_reactants_at_surfaces);
    
    result2.rate_constant_scale = DimensionedExpression(1);
    % result2.rate_constant_scale = result2.rate_constant_scale * parameterReference('N_Avo')^(given_n_reactants - 1);
    result2.rate_constant_scale = result2.rate_constant_scale * parameterReference('N_Avo')^(result2.volume_exponent);
    % result2.rate_constant_scale = result2.rate_constant_scale * given_compartment_volume^(given_n_reactants - 1);
    result2.rate_constant_scale = result2.rate_constant_scale * given_compartment_volume^(result2.volume_exponent);
    result2.rate_constant_scale = result2.rate_constant_scale * (given_compartment_volume / parameterReference('eff_width'))^(result2.membrane_exponent);
    
    result2.rate_constant_units = mcell_reaction_rate_constant_unit_function(given_n_reactants, given_n_reactants_at_surfaces);
end

% warning('CellOrganizer:instance2MCellMDL:assumingUnits', 'Assuming options.resolution is in um');
warning('CellOrganizer:instance2MCellMDL', 'Assuming options.resolution is in um');
% resolution = unit_convert('um', mcell_length_unit, options.resolution.cubic);

net_effective_width = DimensionedExpression(options.output.NET.effective_width, net_length_unit);
mcell_effective_width = net_effective_width;
mcell_effective_width_string = 'eff_width';

use_image_adjacency = options.output.NET.use_image_adjacency;

options.output.remove_mesh_intersections = true;
remove_mesh_intersections = options.output.remove_mesh_intersections;
if remove_mesh_intersections
    min_clearance = 0;
    % min_clearance = options.oobbuffer;
    try
        % meshData = removeMeshIntersections(meshData, min_clearance, 'delete_higher_index');
        meshData = removeMeshIntersections(meshData, min_clearance, 'delete_intersects_framework');
    catch ME
        if any(strcmp(ME.identifier, {'CellOrganizer:removeMeshIntersections:selfIntersection', 'CellOrganizer:removeMeshIntersections:intersection'}))
            warning('CellOrganizer:instance2MCellMDL:meshIntersection', 'removeMeshIntersections threw ''%s''', ME.identifier);
            result = false;
            return;
        else
            rethrow(ME);
        end
    end
end



% Constants

mcell_diffusion_coefficient = DimensionedExpression(NETDefaultDiffusionCoefficient, net_diffusion_coefficient_unit);

mesh_name_map = containers.Map('KeyType', 'char', 'ValueType', 'char');
mesh_name_map('cell') = 'PM';
mesh_name_map('nuc') = 'NM';




% Names of loaded models other than cell and nucleus

function x = isStringNumeric(x)
    x = ~isnan(str2double(x));
end




% Minimum and maximum coordinates over all meshes

minimum_coordinates = inf(1, 3);
maximum_coordinates = -inf(1, 3);
for i = 1:length(meshData)
    single_meshData = meshData(i);
    for j = 1:length(single_meshData.list)
        single_meshData_list_item = single_meshData.list(j);
        switch single_meshData_list_item.type
            case 'triangle mesh'
                mesh_vertices = single_meshData_list_item.mesh.vertices;
                minimum_coordinates = min(min(mesh_vertices, [], 1), minimum_coordinates);
                maximum_coordinates = max(max(mesh_vertices, [], 1), maximum_coordinates);
                
            otherwise
                warning('CellOrganizer:instance2MCellMDL:unrecognizedMeshType', 'Unrecognized meshData type ''%s''', single_meshData_list_item.type);
        end
    end
end



% Following based on documentation and a MCellMDL file saved using CellBlender


% If was_given_mcell_file is true, instance2MCellMDL will assume files matching MCellMDLInputFilenamePattern is complete and only attempt to replace geometry
was_given_net_file = ~isempty(NETfile);
was_given_mcell_file = ~isempty(MCellMDLInputFilenamePattern);
was_given_biochemistry_file = was_given_net_file || was_given_mcell_file;
if was_given_mcell_file && was_given_net_file
    error('Giving a NET file and an MCell MDL file together not supported');
end


file_handle = [];

function oprintf(varargin)
    fprintf(file_handle, varargin{:});
end

function fprintfln(varargin)
    if length(varargin) > 1
        fprintf(varargin{1}, [varargin{2}, '\n'], varargin{3:end});
    else
        fprintf(varargin{1}, '\n');
    end
end

function oprintfln(varargin)
    fprintfln(file_handle, varargin{:});
end

function [result2_in_units, result2_units] = de_to_char_with_correct_units(given_expression, given_correct_units)
    given_expression_in_units = given_expression / given_correct_units;
    result2_in_units = char(given_expression_in_units);
    result2_units = char(given_correct_units);
    if ~given_expression_in_units.isDimensionless()
        error('result2_in_units should be dimensionless but has dimensions ''%s''', char(given_expression_in_units.units_manager()));
    end
end

function result2 = final_expression_evaluation_function(given_expression, given_return_double)
    if nargin < 2
        return_double = true;
    end
    result2 = mathEval(given_expression, final_parameters_names_to_expressions_map, struct('return_double', given_return_double));
end


seed_digits = ceil(log10(MCellMDLNumSimulations+1));
seed_format = ['_seed%0', num2str(seed_digits), 'i'];

compartments_species_counts = containers.Map();
compartments_keys = network_info.compartments.keys();
for i = 1:length(network_info.compartments)
    compartments_species_counts(compartments_keys{i}) = 0;
end
for i = 1:length(network_info.species)
    single_species = network_info.species(i);
    compartment = single_species.chosen_compartment;
    compartments_species_counts(compartment) = compartments_species_counts(compartment) + 1;
end


function [result2_base, result2_base_relative] = savepath_function(given_section)
    result2_base_name = [savepath_name, '.', given_section];
    result2_base_relative = [result2_base_name, savepath_ext];
    result2_base = fullfile(savepath_dir, result2_base_relative);
end

savepaths = containers.Map('KeyType', 'char', 'ValueType', 'char');
savepaths_relative = containers.Map('KeyType', 'char', 'ValueType', 'char');
copypaths = containers.Map('KeyType', 'char', 'ValueType', 'char');
copypaths_relative = containers.Map('KeyType', 'char', 'ValueType', 'char');
sections = {'parameters', 'initialization', 'molecules', 'surface_classes' 'reactions', 'geometry', 'geometry_instantiation', 'main'};
sections_to_write = sections;
% sections_to_copy = {};
if was_given_mcell_file
    sections_to_write = {'geometry'};
    temp = cellfun(@(x)~any(strcmp(x, sections_to_write)), sections);
    temp = sections(temp);
    temp = reshape(temp, [], 1);
    % sections_to_copy = [strrep(MCellMDLInputFilenamePattern, '*', temp), temp];
end
for i = 1:length(sections_to_write)
    section = sections_to_write{i};
    [savepaths(section), savepaths_relative(section)] = savepath_function(section);
end
for i = 1:length(sections)
    section = sections{i};
    [copypaths(section), copypaths_relative(section)] = savepath_function(section);
end




function [given_file_contents] = replace_with_transform_function(given_file_contents, given_pattern, given_transform_function)
    % Replace instances of given_pattern in given_file_contents with the return value of given_transform_function
    [given_file_contents_matches_starts, given_file_contents_matches_ends, given_file_contents_matches_tokens] = regexp(given_file_contents, given_pattern, 'start', 'end', 'tokens');
    for j = length(given_file_contents_matches_tokens):-1:1
        given_file_contents_match_tokens = given_file_contents_matches_tokens{j};
        given_file_contents_match_start = given_file_contents_matches_starts(j);
        given_file_contents_match_end = given_file_contents_matches_ends(j);
        
        given_file_contents_match = given_file_contents(given_file_contents_match_start:given_file_contents_match_end);
        given_file_contents_match_tokens2 = given_transform_function(given_file_contents_match_tokens);
        given_file_contents_match2 = strjoin(given_file_contents_match_tokens2, '');
        given_file_contents = [given_file_contents(1:given_file_contents_match_start-1), given_file_contents_match2, given_file_contents(given_file_contents_match_end+1:end)];
    end
end


mcell_include_pattern = '';
mcell_viz_output_pattern = '';
mcell_react_output_pattern = '';
mcell_partition_pattern = '';
function [given_tokens] = mcell_include_pattern_transform_function(given_tokens)
    [~, given_tokens{2}] = savepath_function(given_tokens{2});
end
partitions_use_input_file_step = true;
function [given_tokens] = mcell_partition_dimension_transform(given_tokens)
    given_dimension = given_tokens{3};
    given_i = axis_to_index(given_dimension);
    given_step = str2double(given_tokens{10});
    if ~partitions_use_input_file_step
        given_step = (maximum_coordinates(given_i) - minimum_coordinates(given_i)) / n_partitions;
    end
    given_tokens{5} = sprintf('PARTITION_%s = [[%16e TO %16e STEP %16e]]', given_dimension, minimum_coordinates(given_i), maximum_coordinates(given_i), given_step);
end


if was_given_mcell_file
    % Replace paths in input files
    any_within_line_regexp = '[^\r\n]*?';
    line_start_regexp = '(?:^|\n)';
    line_end_regexp = '(?:$|\n)';
    pattern_regexp = ['(.*)', filesep, '(.*)\.\*(\.mdl)$'];
    temp = regexp(MCellMDLInputFilenamePattern, pattern_regexp, 'tokens');
    temp = temp{1};
    pattern_path = temp{1};
    pattern_prefix = temp{2};
    pattern_extension = temp{3};
    pattern_files_regexp_internal = [regexptranslate('escape', pattern_prefix), '\.(', any_within_line_regexp, ')', regexptranslate('escape', pattern_extension)];
    pattern_files_regexp = ['^', regexptranslate('escape', pattern_path), filesep, pattern_files_regexp_internal, '$'];
    pattern_files = arrayfun(@(x)[pattern_path, filesep, x.name], dir(MCellMDLInputFilenamePattern), 'UniformOutput', false);
    temp = regexprep(pattern_files, pattern_files_regexp, ['$1']);
    pattern_files_copies = cellfun(@(x)savepath_function(x), temp, 'UniformOutput', false);
    
    mcell_include_pattern = ['(', line_start_regexp, '\s*?INCLUDE_FILE\s*?=\s*?")', pattern_files_regexp_internal, '("\s*?)', line_end_regexp];
    float_pattern = '[0-9Ee.+-]+';
    mcell_partition_pattern = ['(', line_start_regexp, '\s*?)(PARTITION_)([XYZ])(\s*=\s*)', '(\[\s*\[\s*)(', float_pattern, ')(\s*TO\s*)(', float_pattern, ')(\s*STEP\s*)(', float_pattern, ')(\s*\]\s*\])', '("\s*?)', line_end_regexp];
    
    for i = 1:length(pattern_files_copies)
        % Read original file's contents
        file_id = fopen(pattern_files{i}, 'r');
        file_contents_lines = textscan(file_id, '%s', 'Delimiter', {'\r\n', '\n', '\r'});
        file_contents_lines = file_contents_lines{1};
        file_contents = strjoin(file_contents_lines, '\n');
        fclose(file_id);
        
        % Modify include paths
        file_contents = replace_with_transform_function(file_contents, mcell_include_pattern, @mcell_include_pattern_transform_function);
        
        % Modify partitions
        file_contents = replace_with_transform_function(file_contents, mcell_partition_pattern, @mcell_partition_dimension_transform);
        
        % Write modified copy
        file_copy_id = fopen(pattern_files_copies{i}, 'w');
        fwrite(file_copy_id, file_contents);
        fclose(file_copy_id);
    end
end



% Parameters MCell MDL file



section = 'parameters';
if savepaths.isKey(section)

    file_handle = fopen(savepaths(section), 'w');



    oprintfln('/*');
    oprintfln('MCell MDL generated by CellOrganizer');
    oprintfln('*/');

    % Variable definition commands

    oprintfln();
    parameters_names_order = network_info.parameters_names_topological_order;
    for i = 1:length(network_info.parameters)
        parameter_name = parameters_names_order{i};
        parameter = network_info.parameters(parameter_name);
        oprintfln('%s = %s /* %s, %s */', parameter_name, char(parameter.value), char(evaluateExpression(parameter)), char(parameter.units_manager));
    end


    fclose(file_handle);

end



% Initialization MCell MDL file



section = 'initialization';
if savepaths.isKey(section)

    file_handle = fopen(savepaths(section), 'w');



    oprintfln('/*');
    oprintfln('MCell MDL generated by CellOrganizer');
    oprintfln('*/');

    % Initialization commands

    oprintfln();
    oprintfln('ACCURATE_3D_REACTIONS = TRUE');
    if ~isnan(MCellMDLInteractionRadius)
        oprintfln('INTERACTION_RADIUS = %e', MCellMDLInteractionRadius);
    else
        error('isnan(MCellMDLInteractionRadius)')
    end
    oprintfln('NOTIFICATIONS');
    oprintfln('{');
    oprintfln('    ALL_NOTIFICATIONS = ON');
    oprintfln('}');
    oprintfln('WARNINGS');
    oprintfln('{');
    oprintfln('    ALL_WARNINGS = WARNING');
    oprintfln('}');


    fclose(file_handle);

end



% Molecules MCell MDL file



section = 'molecules';
if savepaths.isKey(section)

    file_handle = fopen(savepaths(section), 'w');



    oprintfln('/*');
    oprintfln('MCell MDL generated by CellOrganizer');
    oprintfln('*/');

    % Molecule definition commands

    oprintfln();
    for i = 1:length(network_info.species)
        single_species = network_info.species(i);
        single_species_diffusion_coefficient = DimensionedExpression(single_species.diffusion_coefficient);
        oprintfln('DEFINE_MOLECULE %s', single_species.name);
        oprintfln('{');
        if is_membrane_function(single_species.chosen_compartment)
            command_name = 'DIFFUSION_CONSTANT_2D';
        else
            command_name = 'DIFFUSION_CONSTANT_3D';
        end
        [single_species_diffusion_coefficient_in_units_char, single_species_diffusion_coefficient_units_char] = de_to_char_with_correct_units(single_species_diffusion_coefficient, mcell_diffusion_coefficient_unit);
        oprintfln('    %s = %s /* %s, %s */', command_name, single_species_diffusion_coefficient_in_units_char, char(evaluateExpression(single_species_diffusion_coefficient)), single_species_diffusion_coefficient_units_char);
        oprintfln('}');
    end


    fclose(file_handle);

end



% Surface classes MCell MDL file



section = 'surface_classes';
if savepaths.isKey(section)

    file_handle = fopen(savepaths(section), 'w');



    oprintfln('/*');
    oprintfln('MCell MDL generated by CellOrganizer');
    oprintfln('*/');

    % Geometry definition commands (surface classes)

    for i = 1:length(meshData)
        single_meshData = meshData(i);
        for j = 1:length(single_meshData.list)
            single_meshData_list_item = single_meshData.list(j);
            switch single_meshData_list_item.type
                case 'triangle mesh'
                    surface_class_name = single_meshData_list_item.name;
                    if ~network_info.compartments.isKey(surface_class_name)
                        error('NET file should have extracellular, plasma membrane, cytoplasm, nuclear membrane, and nucleus compartments and membrane and lumen compartments for each protein model. See ''LotkaVolterra.bngl'' for example.');
                    end
                    if network_info.compartments(surface_class_name).spatial_dimensions == 3
                        surface_class_name = network_info.compartments(surface_class_name).outside;
                    end
                    
                    oprintfln();
                    oprintfln('DEFINE_SURFACE_CLASS %s {}', surface_class_name);
                    
                otherwise
                    warning('CellOrganizer:instance2MCellMDL:unrecognizedMeshType', 'Unrecognized meshData type ''%s''', single_meshData_list_item.type);
            end
        end
    end


    fclose(file_handle);

end



% Reactions MCell MDL file



section = 'reactions';
if savepaths.isKey(section)

    file_handle = fopen(savepaths(section), 'w');



    oprintfln('/*');
    oprintfln('MCell MDL generated by CellOrganizer');
    oprintfln('*/');

    % Reaction definition commands

    oprintfln();
    oprintfln('DEFINE_REACTIONS');
    oprintfln('{');
    for i = 1:length(network_info.reactions)
        reaction = network_info.reactions(i);
        reaction_name = reaction.name;
        reaction_rate_constant_net = reaction.rate_constant;
        % reaction_rate_constant_net_units = reaction.rate_constant_units;
        reaction_n_reactants = length(reaction.reactant_indices);
        reaction_compartment = reaction.compartment;
        
        % Check if orientations are necessary
        reactants_in_volume = false;
        reactants_at_surface = false;
        products_in_volume = false;
        products_at_surface = false;
        products_surfaces = {};
        reactants_indices_at_surfaces = [];
        reactants_indices_in_volumes = [];
        products_indices_at_surfaces = [];
        products_indices_in_volumes = [];
        for j = 1:reaction_n_reactants
            if is_membrane_function(network_info.species(reaction.reactant_indices(j)).chosen_compartment)
                reactants_at_surface = true;
                reactants_indices_at_surfaces(end+1) = j;
            else
                reactants_in_volume = true;
                reactants_indices_in_volumes(end+1) = j;
            end
        end
        for j = 1:length(reaction.product_indices)
            if is_membrane_function(network_info.species(reaction.product_indices(j)).chosen_compartment)
                products_at_surface = true;
                products_surfaces{end+1} = network_info.species(reaction.product_indices(j)).chosen_compartment;
                products_indices_at_surfaces(end+1) = j;
            else
                products_in_volume = true;
                products_indices_in_volumes(end+1) = j;
            end
        end
        reaction_in_volume = ~reactants_at_surface && ~products_at_surface;
        volume_reactants_surface_product = ~reactants_at_surface && products_at_surface;
        
        reaction.reactants_indices_at_surfaces = reactants_indices_at_surfaces;
        reaction.reactants_indices_in_volumes = reactants_indices_in_volumes;
        reaction.products_indices_at_surfaces = products_indices_at_surfaces;
        reaction.products_indices_in_volumes = products_indices_in_volumes;
        reaction.reaction_in_volume = reaction_in_volume;
        reaction.volume_reactants_surface_product = volume_reactants_surface_product;
        
        orientation_mark = '';
        if ~reaction_in_volume
            orientation_mark = ';';
        end
        
        reactant_surface_location_suffix = '';
        if volume_reactants_surface_product
            % Assume only one surface compartment
            reactant_surface_location_suffix = [' @ ', products_surfaces{1}];
        end
        
        oprintf('    ');
        for j = 1:reaction_n_reactants
            if j > 1
                oprintf(' + ');
            end
            if ~is_membrane_function(network_info.species(reaction.reactant_indices(j)).chosen_compartment) && ~isempty(reactant_surface_location_suffix)
                current_reactant_surface_location_suffix = [reactant_surface_location_suffix, orientation_mark];
            else
                current_reactant_surface_location_suffix = '';
            end
            oprintf([network_info.species(reaction.reactant_indices(j)).name, orientation_mark, current_reactant_surface_location_suffix]);
        end
        oprintf(' -> ');
        if length(reaction.product_indices) == 0
            oprintf('NULL');
        end
        for j = 1:length(reaction.product_indices)
            if j > 1
                oprintf(' + ');
            end
            oprintf('%s%s', network_info.species(reaction.product_indices(j)).name, orientation_mark);
        end
        
        
        reaction_rate_info = mcell_reaction_rate_info_function(reaction);
        reaction_rate_constant_units = reaction_rate_info.rate_constant_units;
        reaction_rate_constant_scale = reaction_rate_info.rate_constant_scale;
        reaction_rate_units = reaction_rate_info.rate_units;
        reaction_rate_scale = reaction_rate_info.rate_scale;
        
        reaction_rate_constant = reaction_rate_constant_net * reaction_rate_constant_scale;
        
        use_reaction_rate_hack2 = MCellMDLUseReactionRateHack && reaction_in_volume && reaction_n_reactants == 2;
        if use_reaction_rate_hack2
            warning('*** HACK *** Multiplying reaction_rate_constant by 1e6 for presentation');
            reaction_rate_constant = reaction_rate_constant * 1e6;
        end
        
        [reaction_rate_constant_in_units_char, reaction_rate_constant_units_char] = de_to_char_with_correct_units(reaction_rate_constant, reaction_rate_constant_units);
        
        oprintf(' [ %s ]', reaction_rate_constant_in_units_char);
        oprintf(' : %s', reaction_name);
        oprintf(' /* %s, %s */', char(evaluateExpression(reaction_rate_constant)), reaction_rate_constant_units_char);
        if use_reaction_rate_hack2
            oprintf(' /* *** HACK *** Multiplying reaction_rate_constant by 1e6 for presentation */');
        end
        oprintfln();
    end
    oprintfln('}');


    fclose(file_handle);

end



% Geometry MCell MDL file



info_file_handle = fopen(fullfile(savepath_dir, [savepath_name, '.info.json']), 'w');

info = struct();
info.meshData = meshData;
for i = 1:length(info.meshData)
    info.meshData(i).list = rmfield(meshData(i).list, {'mesh', 'img'});
end
fprintfln(info_file_handle, jsonencode(info));
fclose(info_file_handle);



section = 'geometry';
if savepaths.isKey(section)

    file_handle = fopen(savepaths(section), 'w');



    oprintfln('/*');
    oprintfln('MCell MDL generated by CellOrganizer');
    oprintfln('*/');

    % Geometry definition commands (geometrical objects)

    for i = 1:length(meshData)
        single_meshData = meshData(i);
        if options.output.use_individual_meshes && ~strcmp(single_meshData.name, 'frameworkMesh')
            % Combine list into one mesh
            single_meshData_mesh = struct('vertices', zeros(0, 3), 'faces', zeros(0, 3));
            single_meshData_mesh_name = '';
            for j = 1:length(single_meshData.list)
                single_meshData_list_item = single_meshData.list(j);
                
                switch single_meshData_list_item.type
                    case 'triangle mesh'
                        temp = single_meshData_list_item.mesh;
                        single_meshData_mesh_name = single_meshData_list_item.name;
                        single_meshData_mesh.faces = [single_meshData_mesh.faces; temp.faces + size(single_meshData_mesh.vertices, 1)];
                        single_meshData_mesh.vertices = [single_meshData_mesh.vertices; temp.vertices];
                    otherwise
                        warning('CellOrganizer:instance2MCellMDL:unrecognizedMeshType', 'Unrecognized meshData type ''%s''', single_meshData_list_item.type);
                end
            end
            single_meshData.list(1).name = single_meshData_mesh_name;
            single_meshData.list(1).mesh = single_meshData_mesh;
            single_meshData.list = single_meshData.list(1);
        end
        
        for j = 1:length(single_meshData.list)
            single_meshData_list_item = single_meshData.list(j);
            item_name = single_meshData_list_item.name;
            item_type = single_meshData_list_item.type;
            item_mesh = single_meshData_list_item.mesh;
            switch item_type
                case 'triangle mesh'
                    oprintfln();
                    oprintfln('%s POLYGON_LIST', item_name);
                    oprintfln('{');
                    oprintfln('    VERTEX_LIST');
                    oprintfln('    {');
                    for j = 1:size(item_mesh.vertices, 1)
                        vertex = item_mesh.vertices(j, :);
                        oprintfln('        [ %e , %e , %e ]', vertex(1), vertex(2), vertex(3));
                    end
                    oprintfln('    }');
                    oprintfln('    ELEMENT_CONNECTIONS');
                    oprintfln('    {');
                    for j = 1:size(item_mesh.faces, 1)
                        face = item_mesh.faces(j, :) - 1;
                        oprintfln('        [ %i , %i , %i ]', face(1), face(2), face(3));
                    end
                    oprintfln('    }');
                    oprintfln('}');
            end
        end
    end


    fclose(file_handle);

end



% Geometry instantiation MCell MDL file



section = 'geometry_instantiation';
if savepaths.isKey(section)

    file_handle = fopen(savepaths(section), 'w');



    oprintfln('/*');
    oprintfln('MCell MDL generated by CellOrganizer');
    oprintfln('*/');

    % Geometry definition commands (instantiation)

    oprintfln();

    scene_name = 'Scene';
    oprintfln('INSTANTIATE %s OBJECT', scene_name);
    oprintfln('{');

    for i = 1:length(meshData)
        single_meshData = meshData(i);
        for j = 1:length(single_meshData.list)
            single_meshData_list_item = single_meshData.list(j);
            switch single_meshData_list_item.type
                case 'triangle mesh'
                    compartment_object_name = single_meshData_list_item.name;
                    oprintfln();
                    oprintfln('    %s OBJECT %s {}', compartment_object_name, compartment_object_name);
                    
                otherwise
                    warning('CellOrganizer:instance2MCellMDL:unrecognizedMeshType', 'Unrecognized meshData type ''%s''', single_meshData_list_item.type);
            end
        end
    end

    oprintfln();
    for i = 1:length(network_info.species)
        single_species = network_info.species(i);
        single_species_name = single_species.name;
        single_species_release_site_name = [single_species_name, '_release_site'];
        single_species_compartment = single_species.chosen_compartment;
        single_species_compartment = single_species.chosen_compartment;
        single_species_concentration = single_species.concentration;
        
        volume_to_expression_function = @(x)sprintf('%s.%s[ALL]', scene_name, x);
        membrane_to_expression_function = @(x)sprintf('%s.%s[ALL]', scene_name, x);
        
        if is_membrane_function(single_species_compartment)
            orientation_mark = ';';
            shape_expression = membrane_to_expression_function(single_species_compartment);
        else
            orientation_mark = '';
            shape_expression = volume_to_expression_function(network_info.compartments(single_species_compartment).outside);
            single_species_compartment_insides = network_info.compartments(single_species_compartment).insides;
            single_species_compartment_all_insides_insides = {};
            for j = 1:length(single_species_compartment_insides)
                single_species_compartment_all_insides_insides = [single_species_compartment_all_insides_insides, network_info.compartments(single_species_compartment_insides{j}).insides];
            end
            single_species_compartment_all_insides_insides = unique(single_species_compartment_all_insides_insides);
                compartment_surface_class_names = {surface_class_name};
            for j = 1:length(single_species_compartment_all_insides_insides)
                if j == 1
                    shape_expression = [shape_expression, ' - ('];
                else
                    shape_expression = [shape_expression, ' + '];
                end
                shape_expression = [shape_expression, volume_to_expression_function(network_info.compartments(single_species_compartment_all_insides_insides{j}).outside)];
                if j == length(single_species_compartment_all_insides_insides)
                    shape_expression = [shape_expression, ')'];
                end
            end
        end
        
        oprintfln('    %s RELEASE_SITE', single_species_release_site_name);
        oprintfln('    {');
        oprintfln('        SHAPE = %s', shape_expression);
        oprintfln('        MOLECULE = %s%s', single_species_name, orientation_mark);
        
        mcell_single_species_concentration = single_species_concentration;
        
        if is_membrane_function(single_species_compartment)
            mcell_single_species_density = mcell_single_species_concentration * parameterReference('N_Avo') * parameterReference('eff_width');
            
            [mcell_single_species_density_in_units_char, mcell_single_species_density_units_char] = de_to_char_with_correct_units(mcell_single_species_density, mcell_density_unit);
            oprintfln('        DENSITY = %s /* %s, %s */', mcell_single_species_density_in_units_char, char(evaluateExpression(mcell_single_species_density)), mcell_single_species_density_units_char);
        else
            [mcell_single_species_concentration_in_units_char, mcell_single_species_concentration_units_char] = de_to_char_with_correct_units(mcell_single_species_concentration, mcell_concentration_unit);
            oprintfln('        CONCENTRATION = %s /* %s, %s */', mcell_single_species_concentration_in_units_char, char(evaluateExpression(mcell_single_species_concentration)), mcell_single_species_concentration_units_char);
        end
        oprintfln('        RELEASE_PROBABILITY = %s', num2str(1));
        oprintfln('    }');
    end

    oprintfln('}');


    fclose(file_handle);

end



% Seeded simulation MCell MDL file


savepath_viz_output_name = fullfile('viz_data', savepath_name);
savepath_viz_output_relative = [savepath_viz_output_name, savepath_ext];
savepath_viz_output = fullfile(savepath_dir, savepath_viz_output_relative);
savepath_react_output_name = fullfile('react_data', savepath_name);
savepath_react_output_relative = [savepath_react_output_name, savepath_ext];
savepath_react_output = fullfile(savepath_dir, savepath_react_output_relative);


section = 'main';
if savepaths.isKey(section)

    file_handle = fopen(savepath, 'w');


    oprintfln('/*');
    oprintfln('MCell MDL generated by CellOrganizer');
    oprintfln('*/');


    % Initialization commands

    oprintfln();
    oprintfln('TIME_STEP = %e', MCellMDLDefaultTimeStep);
    oprintfln('TIME_STEP_MAX = %e', MCellMDLMaxTimeStep);
    oprintfln('MICROSCOPIC_REVERSIBILITY = ON');
    oprintfln('ITERATIONS = %i', ceil(MCellMDLEndTime / MCellMDLDefaultTimeStep));
    oprintfln('CHECKPOINT_REALTIME = %s', MCellMDLMaxRealTime);


    % Set up partitions

    for i = 1:3
        oprintfln('PARTITION_%s = [ [%e TO %e STEP %e] ]', index_to_axis{i}, minimum_coordinates(i), maximum_coordinates(i), (maximum_coordinates(i) - minimum_coordinates(i)) / n_partitions);
    end


    % Include commands

    oprintfln();
    for i = 1:length(MCellMDLInputFilenamePattern)
        oprintfln('INCLUDE_FILE = "%s"', savepaths_relative(MCellMDLInputFilenamePattern{i}));
    end
    for i = 1:length(sections_to_write)
        oprintfln('INCLUDE_FILE = "%s"', savepaths_relative(sections_to_write{i}));
    end


    % Output specification commands (visualization output)

    oprintfln();
    oprintfln('sprintf(seed,"%%05g",SEED)');


    % Output specification commands (visualization output)

    oprintfln();
    oprintfln('VIZ_OUTPUT');
    oprintfln('{');
    oprintfln('    MODE = CELLBLENDER');
    oprintfln('    FILENAME = "%s/seed_" & seed', savepath_viz_output_name);
    % oprintfln('    FILENAME = "%s"', savepath_viz_output);
    oprintfln('    MOLECULES');
    oprintfln('    {');
    oprintfln('        NAME_LIST { ALL_MOLECULES }');
    oprintfln('        TIME_POINTS { POSITIONS @ [ [ 0 TO %e STEP %e ] ] }', MCellMDLEndTime, MCellMDLOutputTimeStep);
    oprintfln('    }');
    oprintfln('}');


    % Output specification commands (reaction output)

    oprintfln();
    oprintfln('REACTION_DATA_OUTPUT');
    oprintfln('{');
    oprintfln('    STEP = %e', MCellMDLOutputTimeStep);
    for i = 1:length(network_info.species)
        single_species = network_info.species(i);
        single_species_name = single_species.name;
        oprintfln('    {COUNT[%s,WORLD]} => "%s/seed_" & seed & "/%s.World.dat"', single_species_name, savepath_react_output_name, single_species_name);
    end
    for i = 1:length(network_info.reactions)
        reaction = network_info.reactions(i);
        reaction_name = reaction.name;
        oprintfln('    {COUNT[%s,WORLD]} => "%s/seed_" & seed & "/%s.World.dat"', reaction_name, savepath_react_output_name, reaction_name);
    end
    oprintfln('}');


    fclose(file_handle);

end

%{
for i = 1:length(sections_to_copy)
    copyfile(sections_to_copy{i, 1}, copypaths_relative(sections_to_copy{i, 2}));
end
%}



% warning('CellOrganizer:instance2MCellMDL:todo', 'TODO: Double check unit adjustments (um_to_vcml_length_unit, etc.)');
warning('CellOrganizer:instance2MCellMDL', 'TODO: Double check unit adjustments (um_to_vcml_length_unit, etc.)');



result = true;

end

