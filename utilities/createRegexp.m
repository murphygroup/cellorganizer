function pattern = createRegexp(variables_info, options)
%CREATEREGEXP Creates a regular expression pattern from a specification.
%
% Inputs
% ------
% variables_info = struct array containing the following optional fields (default values in parentheses):
%     name ('')
%     type ('any')
%     required (true)
%     pattern ('')
% options        = optional struct array containing the following optional fields (default values in parentheses):
%     separator ('\s')
%     separator_has_variable_length (true)
%     begins_line (false)
%     ends_line (false)
% 
% Outputs
% -------
% pattern = regular expression pattern


% Copyright (C) 2018 Taraz Buck
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.


if nargin < 2
    options = struct();
end
if ~isfield(options, 'separator')
    options.separator = '\s';
end
if ~isfield(options, 'separator_has_variable_length')
    options.separator_has_variable_length = true;
end
if ~isfield(options, 'begins_line')
    options.begins_line = false;
end
if ~isfield(options, 'ends_line')
    options.ends_line = false;
end

if ~isfield(variables_info, 'required')
    variables_info.required = true;
end

if options.begins_line
    pattern = '^';
else
    pattern = '';
end

if options.separator_has_variable_length
    separator_pattern = [options.separator, '+'];
else
    separator_pattern = options.separator;
end

variables_required_indices = find(cell2mat({variables_info.required}));
variables_required_first_index = variables_required_indices(1);
variables_required_last_index = variables_required_indices(end);

for i = 1:length(variables_info)
    variable_info = variables_info(i);
    if ~isfield(variable_info, 'name')
        variable_info.name = '';
    end
    if ~isfield(variable_info, 'type')
        variable_info.type = 'any';
    end
    %{
    if ~isfield(variable_info, 'required')
        variable_info.required = true;
    end
    %}
    if ~isfield(variable_info, 'pattern')
        variable_info.pattern = '';
    end
    
    variable_pattern = '';
    
    variable_pattern = [variable_pattern, '(?:'];
    
    switch variable_info.type
        case 'any'
            type_pattern = ['[^', options.separator, ']'];
        case 'integer'
            type_pattern = '[\-\+]?[0-9]+';
        case 'integer_signless'
            type_pattern = '[0-9]+';
        case 'float'
            type_pattern = '[\-\+]?(?:(?:(?:(?:[0-9]+\.[0-9]*)|(?:[0-9]*\.[0-9]+))(?:[eE][\-\+]?[0-9]+)?)|(?:[iI][nN][fF]))|(?:[nN][aA][nN]))';
        case 'identifier'
            type_pattern = '[a-zA-Z0-9_]+';
        case 'expression'
            type_pattern = '[a-zA-Z0-9_+\-*/\^\(\)\.\s]+';
        case 'expression_spaceless'
            type_pattern = '[a-zA-Z0-9_+\-*/\^\(\)\.]+';
        case 'custom'
            type_pattern = variable_info.pattern;
        otherwise
            error('variable_info.type must be in {''any'', ''integer'', ''float'', ''identifier'', ''expression'', ''custom''}');
    end
    variable_pattern = [variable_pattern, type_pattern];
    
    variable_pattern = [variable_pattern, ')'];
    
    if i < variables_required_last_index
        % Separator must follow
        variable_pattern = ['(?:', variable_pattern, separator_pattern, ')'];
    elseif i == variables_required_last_index
        % Separator already precedes or this is the first token
        variable_pattern = ['(?:', variable_pattern, ')'];
    else
        % Separator must precede
        variable_pattern = ['(?:', separator_pattern, variable_pattern, ')'];
    end
    
    if length(variable_info.name) > 0
        variable_pattern = ['(?<', variable_info.name, '>', variable_pattern, ')'];
    end
    
    if ~variable_info.required
        variable_pattern = [variable_pattern, '?'];
    end
    
    pattern = [pattern, variable_pattern];
end

if options.ends_line
    pattern = [pattern, '$'];
end

    
end
