function network_info = readNetwork(NETfile)
%READNETWORK Parses compartmental BioNetGen-generated NET file.
%
% Inputs
% ------
% NETfile = input filename
%
% Outputs
% -------
% network_info = struct containing parsed NET file data
%
% Notes
% -----
% * NETfile should be the network file generated from a cBNG (compartmental BioNetGen) file by BNG2.pl.


% Authors: Taraz Buck
%
% Copyright (C) 2019 Murphy Lab
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

% 2019-06-11 Taraz Buck - Copied `instance2VCML.m` to `readNetwork.m`.


parameters_names_NET_order = {};
parameters = containers.Map('KeyType', 'char', 'ValueType', 'any');
compartments = containers.Map('KeyType', 'char', 'ValueType', 'any');
molecule_types = struct('index', [], 'type', {}, 'comment', {});
observables = struct('index', [], 'type', {}, 'name', {}, 'observable_expressions', {}, 'comment', {});
species = struct('index', [], 'species_graph', {}, 'count', {}, 'compartments', {}, 'comment', {});
reaction_rules = struct('name', {}, 'reaction', {}, 'direction', {}, 'forward_rate_constant', {}, 'reverse_rate_constant', {}, 'comment', {});
reactions = struct('index', {}, 'reactant_indices', {}, 'product_indices', {}, 'rate_constant', {}, 'unit_conversion', {}, 'comment', {});

if isempty(NETfile) || ~ischar(NETfile)
    % Produce skeleton anyway
    contents = '';
else
    contents = fileread(NETfile);
end
contents_lines = strsplit(contents, {'\n', '\r'});
section = '';

warning('CellOrganizer:readNetwork', 'Check ''reverse_rate_constant''')


bng_patterns = getBNGRegexps();

section_patterns_options = struct('begins_line', true, 'ends_line', true);
reaction_rules_pattern_options = struct('begins_line', true, 'ends_line', true, 'separator', '  ', 'separator_has_variable_length', false);

section_patterns = containers.Map('KeyType', 'char', 'ValueType', 'char');
section_patterns('parameters') = createRegexp(struct('name', {'index', 'name', 'value_expression'}, 'type', {'integer_signless', 'identifier', 'expression_spaceless'}, 'required', num2cell(logical([0, 1, 1]))), section_patterns_options);
section_patterns('compartments') = createRegexp(struct('name', {'name', 'spatial_dimensions', 'size_expression', 'outside'}, 'type', {'identifier', 'integer', 'expression_spaceless', 'identifier'}, 'required', num2cell(logical([1, 1, 0, 0]))), section_patterns_options);
section_patterns('molecule types') = createRegexp(struct('name', {'index', 'type'}, 'type', {'integer_signless', 'custom'}, 'pattern', {'', bng_patterns.net_molecule_type}, 'required', num2cell(logical([0, 1]))), section_patterns_options);
section_patterns('observables') = createRegexp(struct('name', {'index', 'type', 'name', 'observable_expressions'}, 'type', {'integer_signless', 'identifier', 'identifier', 'custom'}, 'pattern', {'', '', '', bng_patterns.net_observables_pattern_list}, 'required', num2cell(logical([0, 1, 1, 1]))), section_patterns_options);
% section_patterns('species') = createRegexp(struct('name', {'index', 'species_graph', 'concentration'}, 'type', {'integer_signless', 'custom', 'expression_spaceless'}, 'pattern', {'', bng_patterns.net_species_pattern, ''}, 'required', num2cell(logical([0, 1, 1]))), section_patterns_options);
% section_patterns('species') = createRegexp(struct('name', {'index', 'species_graph', 'concentration'}, 'type', {'integer_signless', 'custom', 'expression_spaceless'}, 'pattern', {'', bng_patterns.net_pre_species_def, ''}, 'required', num2cell(logical([0, 1, 1]))), section_patterns_options);
section_patterns('species') = createRegexp(struct('name', {'index', 'species_graph', 'count'}, 'type', {'integer_signless', 'custom', 'expression_spaceless'}, 'pattern', {'', bng_patterns.net_pre_species_def, ''}, 'required', num2cell(logical([0, 1, 1]))), section_patterns_options);
section_patterns('reaction rules') = createRegexp(struct('name', {'name', 'reaction', 'rate_constants'}, 'type', {'custom', 'custom', 'custom'}, 'pattern', {bng_patterns.net_reaction_rules_name, bng_patterns.net_reaction_rules_rule, bng_patterns.net_reaction_rules_rate_constant}, 'required', num2cell(logical([0, 1, 1]))), reaction_rules_pattern_options);
section_patterns('reactions') = createRegexp(struct('name', {'index', 'reactant_indices', 'product_indices', 'rate_constant'}, 'type', {'custom', 'custom', 'custom', 'custom'}, 'pattern', {bng_patterns.net_reaction_index, bng_patterns.net_reaction_reactant_indices, bng_patterns.net_reaction_product_indices, bng_patterns.net_reaction_rate_constant}, 'required', num2cell(logical([1, 1, 1, 1]))), section_patterns_options);

% section_patterns('groups') = createRegexp(struct());

species_compartments_pattern = ['@' bng_patterns.string_pattern ''];
reaction_comment_unit_conversion_pattern = ['unit_conversion=(' bng_patterns.net_expression ')'];

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
        warning('CellOrganizer:readNetwork', 'Energy patterns ignored');
        continue;
    elseif strcmp(contents_line, 'end energy patterns')
        section = '';
        continue;
    elseif strcmp(contents_line, 'begin functions')
        section = 'functions';
        warning('CellOrganizer:readNetwork', 'Functions ignored');
        continue;
    elseif strcmp(contents_line, 'end functions')
        section = '';
        continue;
    elseif any(strcmp(contents_line, {'begin species', 'begin seed species'}))
        section = 'species';
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
            info_struct.value_expression = DimensionedExpression(info_struct.value_expression);
            info_struct.comment = contents_line_comment;
            parameters(info_struct.name) = info_struct;
            parameters_names_NET_order{end+1} = info_struct.name;
            
        case 'compartments'
            info_struct = regexp(contents_line, section_patterns('compartments'), 'names');
            info_struct.name = strtrim(info_struct.name);
            info_struct.spatial_dimensions = str2double(info_struct.spatial_dimensions);
            % From Sekar and Faeder 2012:
            % "We considered membrane thickness in dm and surface area in dm2. By simply multiplying these values, we directly get the volume in dm3, which is equivalent to liters."
            info_struct.size_expression = strtrim(info_struct.size_expression);
            info_struct.size_expression = DimensionedExpression(info_struct.size_expression, 'dm3');
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
            molecule_types(end+1) = orderfields(info_struct, molecule_types);
            
        case 'observables'
            % Do nothing
            % warning('CellOrganizer:readNetwork', 'Unfinished');
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
            % info_struct.concentration = strtrim(info_struct.concentration);
            info_struct.count = strtrim(info_struct.count);
            info_struct.count = DimensionedExpression(info_struct.count, 'molecule');
            info_struct.compartments = regexp(info_struct.species_graph, species_compartments_pattern, 'match');
            info_struct.compartments = strrep(info_struct.compartments, '@', '');
            info_struct.compartments = unique(info_struct.compartments);
            % info_struct.concentration = [info_struct.count, '/N_Avo/(', compartments(info_struct.name).size_expression, ')'];
            
            info_struct.comment = contents_line_comment;
            species(end+1) = orderfields(info_struct, species);
            
        case 'reaction rules'
            % Do nothing
            % warning('CellOrganizer:readNetwork', 'Unfinished');
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
            info_struct.unit_conversion = regexp(contents_line_comment, reaction_comment_unit_conversion_pattern, 'tokens');
            if length(info_struct.unit_conversion) > 0
                info_struct.unit_conversion = info_struct.unit_conversion{1}{1};
            else
                info_struct.unit_conversion = '1';
            end
            info_struct.unit_conversion = DimensionedExpression(info_struct.unit_conversion);
            
            info_struct.comment = contents_line_comment;
            
            % From Sekar and Faeder 2012:
            % "Thus, when writing a unimolecular reaction rule in BNGL, one does not need to worry about converting the unimolecular reaction rate constant as long as it is in per second"
            % "When writing a bimolecular reaction rule in BNGL, however, one should convert the bimolecular reaction rate constant to per second"
            % % rate_constant_tokens = strrep(info_struct.rate_constant, 's-1');
            % info_struct.rate_constant = DimensionedExpression(info_struct.rate_constant, 's-1');
            % info_struct.rate_constant = DimensionedExpression(info_struct.rate_constant, 'molecule.s-1');
            % info_struct.rate_constant = DimensionedExpression(info_struct.rate_constant, 1 / (UnitsManager('molecule')^(length(info_struct.reactant_indices) - 1) * UnitsManager('s')));
            info_struct.rate_constant = DimensionedExpression(info_struct.rate_constant, UnitsManager('molecule')^(1 - length(info_struct.reactant_indices)) * UnitsManager('s-1'));
            
            reactions(end+1) = orderfields(info_struct, reactions);
            
        case 'groups'
            % Do nothing
            % warning('CellOrganizer:readNetwork', 'Unfinished');
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

% parameters, compartments, molecule_types, observables, species, reaction_rules, reactions, groups

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

% Add insides to compartments
n_network_info_compartments = length(compartments);
network_info_compartments_keys = compartments.keys;
for i = 1:n_network_info_compartments
    compartment_name = network_info_compartments_keys{i};
    compartment = compartments(compartment_name);
    compartment.insides = {};
    compartments(compartment_name) = compartment;
end

for i = 1:n_network_info_compartments
    compartment_name = network_info_compartments_keys{i};
    compartment = compartments(compartment_name);
    compartment_outside_name = compartment.outside;
    if ~isempty(compartment_outside_name)
        compartment_outside = compartments(compartment_outside_name);
        compartment_outside.insides{end+1} = compartment_name;
        compartments(compartment_outside_name) = compartment_outside;
    end
end

% Add double insides and outsides to compartments

all_compartments_double_insides = containers.Map('KeyType', 'char', 'ValueType', 'any');
all_compartments_double_outsides = containers.Map('KeyType', 'char', 'ValueType', 'any');
for i = 1:length(compartments)
    all_compartments_double_insides(network_info_compartments_keys{i}) = {};
    all_compartments_double_outsides(network_info_compartments_keys{i}) = {};
end

for i = 1:length(compartments)
    compartment = compartments(network_info_compartments_keys{i});
    compartment_name = compartment.name;
    for j = 1:length(compartment.insides)
        compartment_inside = compartment.insides{j};
        all_compartments_double_insides(compartment_name) = [all_compartments_double_insides(compartment_name), compartments(compartment_inside).insides];
    end
    all_compartments_double_insides(compartment_name) = unique(all_compartments_double_insides(compartment_name));
    if ~isempty(compartment.outside)
        all_compartments_double_outsides(compartment_name) = compartments(compartment.outside).outside;
    end
end

for i = 1:n_network_info_compartments
    compartment_name = network_info_compartments_keys{i};
    compartment = compartments(compartment_name);
    compartment.double_outside = all_compartments_double_outsides(compartment_name);
    compartment.double_insides = all_compartments_double_insides(compartment_name);
    compartments(compartment_name) = compartment;
end


% Create return structure

network_info = struct();
network_info.parameters = parameters;
network_info.parameters_names_NET_order = parameters_names_NET_order;
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
network_info.regexp_patterns.reaction_comment_unit_conversion_pattern = reaction_comment_unit_conversion_pattern;

end
