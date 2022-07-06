function [ result ] = instance2VCML( CSGdata, meshData, models, imgs, options, savepath )
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

% 2018-10-02 Taraz Buck - Copied `instance2SBML3_mod.m` to `instance2VCML.m`.


debug = false;
% debug = true;

if ~debug
    warning('off', 'CellOrganizer:instance2VCML');
end

result = false;

if nargin < 3
    % warning('CellOrganizer:instance2VCML:missingRequiredArgument', 'First three arguments are required');
    warning('CellOrganizer:instance2VCML', 'First three arguments are required');
    return;
end

if nargin < 4
    % warning('CellOrganizer:instance2VCML:missingOptionalArgument', 'Argument savepath not given, defaulting to ''./model.vcml''');
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
    % warning('CellOrganizer:instance2VCML:imgsEmpty', 'Argument imgs must contain nonempty 3D arrays');
    warning('CellOrganizer:instance2VCML', 'Argument imgs must contain nonempty 3D arrays');
    return;
end


% Process NET file and combine with geometry
options.output.NET.output_length_unit = 'um';
options.output.NET.output_concentration_unit = 'uM';
options.output.NET.output_membrane_substance_unit = 'molecule';
options.output.NET.output_time_unit = 's';
% options.output.NET.use_image_adjacency = true;
options.output.NET.use_image_adjacency = false;
options.output.NET.translations = options.output.VCML.translations;
options.output.NET.downsampling = options.output.VCML.downsampling;

options.output.temp = struct();
options.output.temp.savepath = savepath;
network_with_geometry_info = readNetworkIntoGeometry(CSGdata, meshData, models, imgs, options);

network_info = network_with_geometry_info.network_info;
geometry_info = network_with_geometry_info.geometry_info;
network_with_geometry_info
network_info
geometry_info

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
get_compartment_pair_string = geometry_info.get_compartment_pair_string;
separate_compartment_images = geometry_info.separate_compartment_images;
combine_compartment_images = geometry_info.combine_compartment_images;
get_compartment_image_boundary_value = geometry_info.get_compartment_image_boundary_value;
crop_compartment_image = geometry_info.crop_compartment_image;
combine_resize_compartment_images = geometry_info.combine_resize_compartment_images;
model_names = geometry_info.model_names;
model_cytonuclearflags = geometry_info.model_cytonuclearflags;
model_names_translated = geometry_info.model_names_translated;
name_map = geometry_info.name_map;
framework_compartment_names = geometry_info.framework_compartment_names;
is_framework_compartment_name_function = geometry_info.is_framework_compartment_name_function;
named_imgs = geometry_info.named_imgs;
% named_imgs_volumes = geometry_info.named_imgs_volumes;
named_imgs_cytonuclearflags = geometry_info.named_imgs_cytonuclearflags;
names = geometry_info.names;
% named_imgs_volumes_keys = geometry_info.named_imgs_volumes_keys;
% named_imgs_volumes_keys_sorted = geometry_info.named_imgs_volumes_keys_sorted;
all_compartments_image_before_downsampling = geometry_info.all_compartments_image_before_downsampling;
all_compartments_image = geometry_info.all_compartments_image;
named_imgs_before_downsampling = geometry_info.named_imgs_before_downsampling;
% resolution_before_downsampling = geometry_info.resolution_before_downsampling;
resolution = geometry_info.resolution;
resolution_single = geometry_info.resolution_single;
% adjacent_pairs_before_downsampling = geometry_info.adjacent_pairs_before_downsampling;
adjacent_pairs = geometry_info.adjacent_pairs;
names_sorted = geometry_info.names_sorted;
n_compartments = geometry_info.n_compartments;
% all_compartment_volumes = geometry_info.all_compartment_volumes;
% all_compartment_exclusive_volumes = geometry_info.all_compartment_exclusive_volumes;
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
% parameter_objects = geometry_info.parameter_objects;
species_index_to_name_function = geometry_info.species_index_to_name_function;
% updateSpeciesNames = geometry_info.updateSpeciesNames;
% expression_evaluation_function = geometry_info.expression_evaluation_function;
evaluateExpression = geometry_info.evaluateExpression;

translateWithDefaultIdentity = geometry_info.translateWithDefaultIdentity;



NETfile = options.output.NET.filename;
VCMLfile = options.output.VCML.input_filename;

VCMLTranslations = options.output.VCML.translations;
VCMLDownsample = options.output.VCML.downsampling;
VCMLAddTranslocationIntermediates = options.output.VCML.add_translocation_intermediates;
% VCMLObjectsAlwaysPresent = options.output.VCML.objects_always_present;

VCMLNumSimulations = 1;
% VCMLNumSimulations = options.output.VCML.num_simulations;
VCMLDeleteInputSimulations = options.output.VCML.delete_input_simulations;
VCMLEndTime = options.output.VCML.end_time;
VCMLDefaultTimeStep = options.output.VCML.default_time_step;
VCMLMinTimeStep = options.output.VCML.min_time_step;
VCMLMaxTimeStep = options.output.VCML.max_time_step;
VCMLOutputTimeStep = options.output.VCML.output_time_step;
VCMLAbsoluteTolerance = options.output.VCML.absolute_tolerance;
VCMLRelativeTolerance = options.output.VCML.relative_tolerance;
VCMLDefaultDiffusionCoefficient = options.output.VCML.default_diffusion_coefficient;

net_length_unit = options.output.NET.units.length;
net_area_unit = [net_length_unit,'2'];
net_volume_unit = [net_length_unit,'3'];
net_time_unit = options.output.NET.units.time;
net_count_unit = 'molecules';
net_concentration_unit = options.output.NET.units.concentration;
% warning('CellOrganizer:instance2VCML:assumingUnits', 'Assuming NET file concentrations are extensive');
warning('CellOrganizer:instance2VCML', 'Assuming NET file concentrations are extensive');

% avogadro_constant = 6.022140857e23;

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

%{
% warning('CellOrganizer:instance2VCML:assumingUnits', 'Assuming options.resolution is in um');
warning('CellOrganizer:instance2VCML', 'Assuming options.resolution is in um');
resolution = unit_convert('um', vcml_length_unit, options.resolution.cubic);
%}

SI_effective_width = options.output.NET.effective_width;
vcml_effective_width = unit_convert('m', vcml_length_unit, SI_effective_width);



% Constants

compartment_volume_prefix = 'vol_';
compartment_area_prefix = 'sa_';
compartment_diffusion_coefficient_prefix = 'dc_';

% kinetics_type = 'GeneralKinetics';
kinetics_type = 'MassAction';

% avogadro_constant_value_expression = '6.022140857e23';
% avogadro_constant_units = 'molecules.mol-1';

vcml_diffusion_coefficient = unit_convert('m2.s-1', vcml_diffusion_coefficient_unit, VCMLDefaultDiffusionCoefficient);




% Process BioNetGen NET file

% network_info = readNetwork(NETfile);


% Helper functions

function x = isStringNumeric(x)
    x = ~isnan(str2double(x));
end



% Create a single image containing all compartments. Assumes no overlapping.

% Collect names and properties of compartments



% Create compartment adjacency matrices for paths of length 1 to adjacent_pairs_max_degree


% Find all objects in each compartment

%{
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
    % warning('CellOrganizer:instance2VCML:assumingOneCell', 'Assuming there is only one cell');
    warning('CellOrganizer:instance2VCML', 'Assuming there is only one cell');
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
%}



% warning('CellOrganizer:instance2VCML:todo', 'Some inconsistencies! Voxel not 6-connected then determined to be 6-adjacent');
warning('CellOrganizer:instance2VCML', 'Some inconsistencies! Voxel not 6-connected then determined to be 6-adjacent');


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
    
    % warning('CellOrganizer:instance2VCML:todo', 'Make values compatible with unit system or make assumptions!');
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

all_compartment_VCML_data = all_compartment_data;
all_object_VCML_data = all_object_data;
all_object_membrane_VCML_data = all_object_membrane_data;
all_membrane_VCML_data = all_membrane_data;


parameters_names_to_expressions_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
function result2 = expression_evaluation_function(given_expression, given_return_double)
    if nargin < 2
        return_double = true;
    end
    result2 = mathEval(given_expression, parameters_names_to_expressions_map, struct('return_double', given_return_double));
end

%{
parameter_objects = struct('index', {}, 'name', {}, 'value', {}, 'is_numeric', {}, 'comment', {}, 'units', {});
if was_given_net_file
    updateSpeciesNames();

    % Add or replace volume and surface area parameters
    network_info_parameters_keys = network_info.parameters.keys;
    % warning('CellOrganizer:instance2VCML:assumingUnits', 'Assuming all unknown parameters have no units');
    warning('CellOrganizer:instance2VCML', 'Assuming all unknown parameters have no units');
    for i = 1:length(network_info.parameters)
        parameter_name = network_info_parameters_keys{i};
        parameter = network_info.parameters(parameter_name);
        parameter.units = '1';
        network_info.parameters(parameter_name) = parameter;
    end

    % Set reaction rate constant units
    network_info_parameters_keys = network_info.parameters.keys;
    % warning('CellOrganizer:instance2VCML:assumingUnits', 'Assuming k_*, kp_*, km_* parameters are reaction rate constants');
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
    % warning('CellOrganizer:instance2VCML:assumingUnits', 'Assuming *_init_count parameters have units ''%s''', net_count_unit);
    warning('CellOrganizer:instance2VCML', 'Assuming *_init_count parameters have units ''%s''', net_count_unit);
    % warning('CellOrganizer:instance2VCML:assumingUnits', 'Assuming *_init_conc parameters have units ''%s''', net_concentration_unit);
    warning('CellOrganizer:instance2VCML', 'Assuming *_init_conc parameters have units ''%s''', net_concentration_unit);
    % warning('CellOrganizer:instance2VCML:assumingValue', 'Assuming N_A parameter is Avogadro constant, has value ''%s'', and has units ''%s''', avogadro_constant_value_expression, avogadro_constant_units);
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
%}



% Set up VCML structure


DiagramNodesMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
DiagramNodesIndices = containers.Map('KeyType', 'char', 'ValueType', 'double');

% warning('CellOrganizer:instance2VCML:todo', 'Ignoring CSGdata.primitiveOnly.');
warning('CellOrganizer:instance2VCML', 'Ignoring CSGdata.primitiveOnly.');
% warning('CellOrganizer:instance2VCML:todo', 'Ordinals not used here, connectivity inferred using adjacency in image.');
warning('CellOrganizer:instance2VCML', 'Ordinals not used here, connectivity inferred using adjacency in image.');


% warning('CellOrganizer:instance2VCML:todo', 'consensus_compartment: Check if this follows rules in Sekar and Faeder 2012');
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


%{
% Decide on the reaction rate and units for each reaction
% warning('CellOrganizer:instance2VCML:assumingUnits', 'network_info.reactions.rate_constant_units values not known to be correct for BNG and/or VCell');
warning('CellOrganizer:instance2VCML', 'network_info.reactions.rate_constant_units values not known to be correct for BNG and/or VCell');
% warning('CellOrganizer:instance2VCML:assumingUnits', 'network_info.reactions.rate_constant values are not processed as microscopic');
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
%}

%{
% Split BNGL translocations because VCell does not accept reactions between non-adjacent compartments
reactions_to_add = network_info.reactions(1:0);
reactions_to_remove = [];
% warning('CellOrganizer:instance2VCML:assumption', 'Splitting BNGL translocations with assumptions about diffusion coefficients and reaction rate, see code');
warning('CellOrganizer:instance2VCML', 'Splitting BNGL translocations with assumptions about diffusion coefficients and reaction rate, see code');
% warning('CellOrganizer:instance2VCML:assumption', 'Permitting BNGL translocations that are invalid according to Sekar and Faeder 2012, "Rule-Based Modeling of Signal Transduction: A Primer" but used in the published example "journal.pcbi.1004611.s003.bngl" with assumptions about diffusion coefficients and reaction rate, see code');
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
%}


if ~was_given_vcml_file
    % Was commented because not known to be necessary. `ModelParameters` does not appear in a VCML file built entirely in VCell itself. Guessing it is generated by SBML import.
    ModelParametersNode = docNode.createElement('ModelParameters');
    parameters_names_order = network_info.parameters_names_topological_order;
    parameters_should_write_to_parameters = network_info.parameters_should_write_to_parameters;
    set_should_write_to_parameters = network_info.set_should_write_to_parameters;
    get_should_write_to_parameters = network_info.get_should_write_to_parameters;
    for i = 1:length(network_info.parameters)
        parameter_name = parameters_names_order{i};
        if ~get_should_write_to_parameters(parameter_name, parameters_should_write_to_parameters)
            continue
        end
        parameter = network_info.parameters(parameter_name);
        parameter_units = parameter.units_manager.toString(false, 'VCML');
        ParameterNode = docNode.createElement('Parameter');
        ParameterNode.setAttribute('Name',parameter_name);
        ParameterNode.setAttribute('Role','user defined');
        ParameterNode.setAttribute('Unit', parameter.units_manager.toString(false, 'VCML'));
        ParameterNode.setTextContent(char(parameter.value));
        ModelParametersNode.appendChild(ParameterNode);
    end
    ModelNode.appendChild(ModelParametersNode);

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
        MembraneNode.setAttribute('MembraneVoltage',object.MembraneVoltageName);
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



    % warning('CellOrganizer:instance2VCML:todo', 'Rate constant units converted improperly');
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
    % warning('CellOrganizer:instance2VCML:assumingUnits', 'ParameterNode with role ''user defined'' attribute Unit not computed, assuming ??? for 1 parameter');
    warning('CellOrganizer:instance2VCML', 'ParameterNode with role ''user defined'' attribute Unit not computed, assuming ??? for 1 parameter');
    [network_info.species.rate] = species_rates{:};
    [network_info.reactions.mass_action_rate] = reaction_mass_action_rates{:};
end



for i = 1:length(all_compartment_data)
    object = all_compartment_data(i);
    name = object.name;
    
    DiagramNode = createOrGetChild(ModelNode, 'Diagram', struct('Name', object.DiagramNodeName, 'Structure', object.FeatureName));
    
    DiagramNodesMap(object.FeatureName) = DiagramNode;
    DiagramNodesIndices(object.FeatureName) = length(DiagramNodesMap);
end

for k = 1:length(all_membrane_data)
    object = all_membrane_data(k);
    i = object.i;
    j = object.j;
    
    start_obj = all_compartment_data(i);
    link_obj = all_compartment_data(j);
    
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
    % warning('CellOrganizer:instance2VCML:todo', 'Not checking ModelUnitSystem');
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
    % warning('CellOrganizer:instance2VCML:todo', 'Unit conversions unfinished');
    warning('CellOrganizer:instance2VCML', 'Unit conversions unfinished');
    % warning('CellOrganizer:instance2VCML:todo', 'Check VolumeSubstanceUnit (''%s'')',vcml_volume_substance_unit);
    warning('CellOrganizer:instance2VCML', 'Check VolumeSubstanceUnit (''%s'')',vcml_volume_substance_unit);
    for i = 1:size(ModelUnitSystemNode_attributes, 1)
        ModelUnitSystemNode.setAttribute(ModelUnitSystemNode_attributes{i, 1}, ModelUnitSystemNode_attributes{i, 2});
    end
    ModelNode.appendChild(ModelUnitSystemNode);
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
            
            GeometryNode2 = getFirstChildByTagName(SimulationSpecNode2, 'Geometry');
            if ~isempty(GeometryNode2)
                GeometryNode2Name = getAttributeValue(GeometryNode2, 'Name', 'str');
                GeometryNode2Name = [options.prefix, '_', GeometryNode2Name];
                GeometryNode2.setAttribute('Name', GeometryNode2Name)
            end
            
            % `GeometryContext` `Name` attribute is not in earlier VCML output that is importable into VCell
            %{
            GeometryContextNode2 = getFirstChildByTagName(SimulationSpecNode2, 'GeometryContext');
            if ~isempty(GeometryContextNode2)
                GeometryContextNode2Name = getAttributeValue(GeometryContextNode2, 'Name', 'str');
                GeometryContextNode2Name = [options.prefix, '_', GeometryContextNode2Name];
                GeometryContextNode2.setAttribute('Name', GeometryContextNode2Name)
            end
            %}
        end
    end
    SimulationSpecNodeName = char(SimulationSpecNode.getAttribute('Name'));
    SimulationNodeName = [SimulationSpecNodeName, '_Simulation_', num2str(i)];
    

    if was_given_vcml_file
        % Translate function domain names
        MathDescriptionNode = getChild(SimulationSpecNode, 'MathDescription');
        if ~isempty(MathDescriptionNode)
            FunctionNodes = MathDescriptionNode.getElementsByTagName('Function');
            for j = 1:FunctionNodes.getLength()
                FunctionNode = FunctionNodes.item(j-1);
                
                % Domain attribute
                FunctionNode_domain = getAttributeValue(FunctionNode, 'Domain', 'char');
                if ~isempty(FunctionNode_domain)
                    FunctionNode_domain2 = FunctionNode_domain;
                    if strends(FunctionNode_domain, '_membrane')
                        FunctionNode_domain2 = FunctionNode_domain2(1:end-length('_membrane'));
                    end
                    FunctionNode_domain2 = translateWithDefaultIdentity(name_map, FunctionNode_domain2);
                    if strends(FunctionNode_domain, '_membrane')
                        FunctionNode_domain2 = [FunctionNode_domain2, '_membrane'];
                    end
                    FunctionNode.setAttribute('Domain', FunctionNode_domain2);
                end
                
                % Text content
                FunctionNode_text = getTextContentValue(FunctionNode, 'char');
                FunctionNode_text_tokens = regexp(FunctionNode_text, '(vcRegionVolume|vcRegionArea)\(''([^'']+)''\)', 'tokens');
                for k = 1:length(FunctionNode_text_tokens)
                    temp = FunctionNode_text_tokens{k, :};
                    FunctionNode_text = strrep(FunctionNode_text, [temp{1}, '(''', temp{2}, ''')'], [temp{1}, '(''', translateWithDefaultIdentity(name_map, temp{2}, {'_membrane'}), ''')']);
                end
                FunctionNode.setTextContent(FunctionNode_text);
            end
        end
    end
    
    
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
    SimulationSpecNode.appendChild(GeometryNode);
        
    extent_ij = size(all_compartments_image);
    extent_ij = extent_ij .* resolution;
    extent_ij_given_vcml = convert_from_vcml_to_given_vcml(extent_ij, 'LengthUnit');
    
    ExtentNode = createOrGetChild(GeometryNode, 'Extent');
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
    for j = 1:length(all_compartment_data)
        object = all_compartment_data(j);
        PixelClassNode = docNode.createElement('PixelClass');
        PixelClassNode.setAttribute('Name',object.PixelClassNodeName);
        PixelClassNode.setAttribute('ImagePixelValue',num2str(j));
        ImageNode.appendChild(PixelClassNode);
    end
    
    removeChildrenByTagName(GeometryNode, 'SubVolume');
    for j = 1:length(all_compartment_data)
        object = all_compartment_data(j);
        SubVolumeNode = docNode.createElement('SubVolume');
        SubVolumeNode.setAttribute('Name',object.SubVolumeNodeName);
        SubVolumeNode.setAttribute('Handle',num2str(j-1));
        SubVolumeNode.setAttribute('Type','Image');
        SubVolumeNode.setAttribute('ImagePixelValue',num2str(j));
        GeometryNode.appendChild(SubVolumeNode);
    end
    
    removeChildrenByTagName(GeometryNode, 'SurfaceClass');
    for j = 1:length(all_membrane_data)
        object = all_membrane_data(j);
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
        object = all_object_data(j);
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
        
    for j = 1:length(all_object_membrane_data)
        object = all_object_membrane_data(j);
        ij = object.ij;
        start_obj = all_object_data(ij(1));
        link_obj = all_object_data(ij(2));
        object_area = object.object_area;
        object_area_given_vcml = convert_from_vcml_to_given_vcml(object_area, 'AreaUnit');
        
        MembraneRegionNode = docNode.createElement('MembraneRegion');
        MembraneRegionNode.setAttribute('Name',object.MembraneRegionNodeName);
        MembraneRegionNode.setAttribute('VolumeRegion1',start_obj.VolumeRegionNodeName);
        MembraneRegionNode.setAttribute('VolumeRegion2',link_obj.VolumeRegionNodeName);
        MembraneRegionNode.setAttribute('Size',num2str(object_area_given_vcml));
        MembraneRegionNode.setAttribute('Unit',vcml_area_unit);
        SurfaceDescriptionNode.appendChild(MembraneRegionNode);
    end
    
    
    % GeometryContext
    
    removeChildrenByTagName(SimulationSpecNode, 'GeometryContext');
    GeometryContextNode = docNode.createElement('GeometryContext');
    SimulationSpecNode.appendChild(GeometryContextNode);
    
    for j = 1:length(all_compartment_data)
        object = all_compartment_data(j);
        object_volume = object.object_volume;
        object_volume_given_vcml = convert_from_vcml_to_given_vcml(object_volume, 'VolumeUnit');
        
        FeatureMappingNode = createOrGetChild(GeometryContextNode, 'FeatureMapping', struct('Feature', object.FeatureName));
        FeatureMappingNode.setAttribute('Feature',object.FeatureName);
        FeatureMappingNode.setAttribute('GeometryClass',object.SubVolumeNodeName);
        FeatureMappingNode.setAttribute('SubVolume',object.SubVolumeNodeName);
        FeatureMappingNode.setAttribute('Size',num2str(object_volume_given_vcml));
        FeatureMappingNode.setAttribute('VolumePerUnitVolume',num2str(1));
    end
    
    for j = 1:length(all_membrane_data)
        object = all_membrane_data(j);
        object_area = object.area;
        object_area_given_vcml = convert_from_vcml_to_given_vcml(object_area, 'AreaUnit');
        
        MembraneMappingNode = createOrGetChild(GeometryContextNode, 'MembraneMapping', struct('Membrane', object.MembraneNodeName));
        MembraneMappingNode.setAttribute('GeometryClass',object.SurfaceClassNodeName);
        MembraneMappingNode.setAttribute('Size',num2str(object_area_given_vcml));
        MembraneMappingNode.setAttribute('AreaPerUnitArea',num2str(1));
        MembraneMappingNode.setAttribute('CalculateVoltage','false');
        MembraneMappingNode.setAttribute('SpecificCapacitance',num2str(1));
        MembraneMappingNode.setAttribute('InitialVoltage',num2str(0));
    end
    
    
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
    
        for k = 1:length(all_membrane_data)
            object = all_membrane_data(k);
            object_area = object.area;
            
            MembraneMappingNode = createOrGetChild(GeometryContextNode2, 'MembraneMapping', struct('Membrane', object.MembraneNodeName));
            MembraneMappingNode.setAttribute('Membrane',object.MembraneNodeName);
            MembraneMappingNode.setAttribute('GeometryClass',object.SurfaceClassNodeName);
            MembraneMappingNode.setAttribute('Size',num2str(object.area));
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
        % DiffusionNode.setTextContent(num2str(object_diffusion_coefficient));
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
    parameters = network_info.parameters;
    parameters_keys = parameters.keys();
    for j = 1:length(parameters)
        parameter_name = parameters_keys{j};
        parameter_object = parameters(parameter_name);
        ParameterNodeName = parameter_name;
        ParameterNodeValue = parameter_object.value;
        parameter_value_is_numeric = ~parameter_object.hasVariables();
        
        if parameter_value_is_numeric
            % Number
            ConstantNode = createOrGetChild(MathDescriptionNode, 'Constant', struct('Name', ParameterNodeName));
            ConstantNode.setTextContent(ParameterNodeValue);
            constant_nodes_to_add{end+1} = ConstantNode;
        else
            % Expression
            FunctionNode = createOrGetChild(MathDescriptionNode, 'Function', struct('Name', ParameterNodeName));
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
    for j = 1:length(all_compartment_data)
        object = all_compartment_data(j);
        
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
    
    
    % warning('CellOrganizer:instance2VCML:todo', 'Unfinished');
    warning('CellOrganizer:instance2VCML', 'Unfinished');
        
    MembraneSubDomainNodeMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    % Keep track of which MembraneSubDomain nodes were in the input file
    temp = MathDescriptionNode.getElementsByTagName('MembraneSubDomain');
    input_MembraneSubDomainNodes = {};
    for j = 1:temp.getLength()
        input_MembraneSubDomainNodes{end+1} = temp.item(j-1);
    end
    
    for j = 1:length(all_membrane_data)
        object = all_membrane_data(j);
        
        MembraneSubDomainNode = createOrGetChild(MathDescriptionNode, 'MembraneSubDomain', struct('InsideCompartment', object.InsideCompartment, 'OutsideCompartment', object.OutsideCompartment));
        MembraneSubDomainNode2 = createOrGetChild(MathDescriptionNode, 'MembraneSubDomain', struct('InsideCompartment', object.OutsideCompartment, 'OutsideCompartment', object.InsideCompartment));
        if ~was_given_biochemistry_file || inCellArray(input_MembraneSubDomainNodes, MembraneSubDomainNode)
            MathDescriptionNode.removeChild(MembraneSubDomainNode2);
        else
            MathDescriptionNode.removeChild(MembraneSubDomainNode);
            MembraneSubDomainNode = MembraneSubDomainNode2;
        end
        
        MembraneSubDomainNode.setAttribute('Name',object.MembraneSubDomainNodeName);
        MembraneSubDomainNodeMap(object.MembraneNodeName) = MembraneSubDomainNode;
        
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
    
    
    % warning('CellOrganizer:instance2VCML:assumingValue', 'JumpCondition always set to zero');
    warning('CellOrganizer:instance2VCML', 'JumpCondition always set to zero');
    % warning('CellOrganizer:instance2VCML:todo', 'Unfinished');
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
            for k = 1:length(all_membrane_data)
                object = all_membrane_data(k);
                if strcmp(object.InsideCompartment, single_species_compartment) || strcmp(object.OutsideCompartment, single_species_compartment)
                    JumpConditionNode = docNode.createElement('JumpCondition');
                    JumpConditionNode.setAttribute('Name',single_species.name);
                    MembraneSubDomainNodeMap(object.MembraneNodeName).appendChild(JumpConditionNode);
                    
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
    
    if was_given_vcml_file && ~VCMLDeleteInputSimulations
        % Copy geometry to other simulations in input file
        GeometryNode_children_to_copy = {
            'Extent',
            'Origin',
            'Image',
            'SubVolume',
            'SurfaceClass',
            'SurfaceDescription',
            };
        GeometryContextNode_children_to_copy = {
            'FeatureMapping',
            'MembraneMapping',
            };
        for j = 1:SimulationSpecNodes.getLength()
            SimulationSpecNode2 = SimulationSpecNodes.item(j-1);
            GeometryNode2 = getFirstChildByTagName(SimulationSpecNode2, 'Geometry');
            GeometryNode2_dimension = getAttributeValue(GeometryNode2, 'Dimension', 'int');
            GeometryContextNode2 = getFirstChildByTagName(SimulationSpecNode2, 'GeometryContext');
            SimulationNode2 = getFirstChildByTagName(SimulationSpecNode2, 'Simulation');
            if GeometryNode2_dimension == 3
                if SimulationSpecNode2 ~= SimulationSpecNode
                    for k = 1:length(GeometryNode_children_to_copy)
                        child_tag_name = GeometryNode_children_to_copy{k};
                        child_node = getFirstChildByTagName(GeometryNode, child_tag_name).copy();
                        removeChildrenByTagName(GeometryNode2, child_tag_name);
                        GeometryNode2.appendChild(child_node);
                    end
                    for k = 1:length(GeometryContextNode_children_to_copy)
                        child_tag_name = GeometryContextNode_children_to_copy{k};
                        child_node = getFirstChildByTagName(GeometryContextNode, child_tag_name).copy();
                        removeChildrenByTagName(GeometryContextNode2, child_tag_name);
                        GeometryContextNode2.appendChild(child_node);
                    end
                end
                
                MeshSpecificationNode2 = getFirstChildByTagName(SimulationNode2, 'MeshSpecification');
                SizeNode2 = getFirstChildByTagName(MeshSpecificationNode2, 'Size');
                SizeNode2.setAttribute('X',num2str(size(all_compartments_image,2)));
                SizeNode2.setAttribute('Y',num2str(size(all_compartments_image,1)));
                SizeNode2.setAttribute('Z',num2str(size(all_compartments_image,3)));
            end
        end
    else
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
        ErrorToleranceNode.setAttribute('Absolut',num2str(VCMLAbsoluteTolerance));
        ErrorToleranceNode.setAttribute('Relative',num2str(VCMLRelativeTolerance));
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
    end
    
    
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
    
    for j = 1:length(all_object_membrane_data)
        VCML_data = all_object_membrane_data(j);
        
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

