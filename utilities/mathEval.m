function value = mathEval(expression, variables_map, options)
%MATHEVAL Evaluates arithmetic expression string without EVAL.
% VALUE = MATHEVAL(EXPRESSION) returns the double value of the evaluated string EXPRESSION or throws an error indicating failure.
%
% Inputs
% ------
% expression    = string containing an arithmetic expression
% variables_map = optional containers.Map from identifier strings to expression strings
% options       = optional struct containing the following optional fields (default values in parentheses) in addition to fields used by mathTokenize:
%     return_type ('double_or_char')
%     verbose (false)
% 
% Outputs
% -------
% value = double value of expression

% Tests
% -----
% assert(mathEval('2 * 3') == 6)
% assert(mathEval('a * b', containers.Map({'a', 'b'}, {'1 + 1', '1 + 1 + 1'})) == 6)
% assert(mathEval('a * b', containers.Map({'a', 'b'}, {'1 + 1', '1 + 1 + 1'}), struct('return_type', 'DimensionedExpression')) == 6)
% assert(mathEval('a * b', containers.Map({'a', 'b'}, {DimensionedExpression(2, 'm'), DimensionedExpression(3, 's-1')}), struct('return_type', 'DimensionedExpression')) == DimensionedExpression(6, 'm.s-1'))
% assert(mathEval('a * b', containers.Map({'a', 'b', 'c', 'd'}, {DimensionedExpression('c + d', 'm'), DimensionedExpression(3, 's-1'), DimensionedExpression(1, 'm'), DimensionedExpression(1, 'm')}), struct('return_type', 'DimensionedExpression')) == DimensionedExpression(6, 'm.s-1'))


% Copyright (C) 2019 Taraz Buck
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


% 2019-01-01 Taraz Buck: Created.


if nargin < 2 || isempty(variables_map)
    variables_map = containers.Map('KeyType', 'char', 'ValueType', 'char');
elseif ~isa(variables_map, 'containers.Map')
    error('CellOrganizer:mathEval', 'variables_map must be containers.Map or empty');
end

if nargin < 3
    options = struct();
end

% Process options
default_options = struct();
default_options.return_type = 'double_or_char';
default_options.verbose = false;

if ~exist('options', 'var')
    options = default_options; 
else
    options = process_options_structure(default_options, options);
end

% Postprocess options

if isa(expression, 'DimensionedExpression')
    original_expression = expression;
    expression = expression.value;
end

%{
if strcmp(options.return_type, 'DimensionedExpression')
    variables_map2 = containers.Map('KeyType', 'char', 'ValueType', 'any');
    variables_map_keys = variables_map.keys();
    for i = 1:length(variables_map)
        variables_map_key = variables_map_keys{i};
        variable = variables_map(variables_map_key);
        if ischar(variable)
            if isStringNumeric(variable)
                variable = str2double(variable);
            else
                variable = DimensionedExpression(variable);
            end
        end
        variables_map2(variables_map_key) = variable;
    end
    variables_map = variables_map2;
    if options.verbose
        for i = 1:length(variables_map)
            variables_map_key = variables_map_keys{i};
            variable = variables_map(variables_map_key);
            fprintf('*** mathEval\n');
            fprintf('***     %s = %s\n', variables_map_key, char(variable));
        end
    end
    for i = 1:length(variables_map2)
        variables_map_key = variables_map_keys{i};
        variable = variables_map2(variables_map_key);
        if isa(variable, 'DimensionedExpression')
            variable_variables = variable.getVariables();
            can_evaluate = true;
            for j = 1:length(variable_variables)
                if ~variables_map2.isKey(variable_variables{j})
                    can_evaluate = false;
                    break;
                end
            end
            can_evaluate = can_evaluate && ~isStringNumeric(variable.value);
            if can_evaluate
                variables_map2.remove(variables_map_key);
                variable = mathEval(variable, variables_map2, options);
                variables_map2(variables_map_key) = variable;
            end
        end
    end
    variables_map = variables_map2;
end
%}

function x = isStringNumeric(x)
    x = ~isnan(str2double(x));
end

function x = valueToTokenValue(x)
    if options.verbose
        fprintf('    valueToTokenValue\n')
        fprintf('    class(x) = %s\n', class(x));
        x
    end
    switch options.return_type
        case {'double_or_char', 'double', 'char'}
            if ischar(x)
                if isStringNumeric(x)
                    x = str2double(x);
                else
                    if variables_map.isKey(x)
                        x = variables_map(x);
                        %{
                        x_units = x.units_manager;
                        x = mathEval(x, variables_map, setfield(options, 'return_type', 'DimensionedExpression'));
                        x = DimensionedExpression(x);
                        % x = DimensionedExpression(x, x_units);
                        %}
                    end
                end
            end
        case 'DimensionedExpression'
            if variables_map.isKey(x)
                x = variables_map(x);
                %{
                x_units = x.units_manager;
                x = mathEval(x, variables_map, setfield(options, 'return_type', 'DimensionedExpression'));
                x = DimensionedExpression(x);
                % x = DimensionedExpression(x, x_units);
                %}
            end
            x = DimensionedExpression(x);
        %{
        case 'char'
            x = x;
        %}
    end
    if options.verbose
        fprintf('    class(x) = %s\n', class(x));
    end
end

function x = valueToReturnValue(x)
    switch options.return_type
        case {'double_or_char', 'double', 'char'}
            if isa(x, 'DimensionedExpression')
                if ~x.isDimensionless()
                    error('Return value is not dimensionless');
                end
                x = x.value;
            end
    end
    switch options.return_type
        case {'double_or_char', 'double'}
            if ischar(x) && isStringNumeric(x)
                x = str2double(x);
            end
        case 'DimensionedExpression'
            if variables_map.isKey(x)
                x = variables_map(x);
                %{
                x_units = x.units_manager;
                x = mathEval(x, variables_map, setfield(options, 'return_type', 'DimensionedExpression'));
                x = DimensionedExpression(x);
                % x = DimensionedExpression(x, x_units);
                %}
            end
            x = DimensionedExpression(x);
        case 'char'
            if isnumeric(x)
                x = double2str(x);
            end
    end
end

function x = valueToString(x)
    switch class(x)
        case 'double'
            x = double2str(x);
        case 'DimensionedExpression'
            x = char(x);
        case 'char'
            x = x;
    end
end


function expressionTokensDisplay(given_expression_tokens)
    given_expression_tokens2 = given_expression_tokens;
    for i = 1:length(given_expression_tokens2)
        if isa(given_expression_tokens2(i).token, 'DimensionedExpression')
            given_expression_tokens2(i).token = char(given_expression_tokens2(i).token);
        end
    end
    given_expression_tokens_cell = [fieldnames(given_expression_tokens2)'; squeeze(struct2cell(given_expression_tokens2))'];
    fprintf('    expressionTokensDisplay:\n');
    disp(given_expression_tokens_cell);
end


% Tokenize

expression_tokens = mathTokenize(expression, variables_map, options);
if options.verbose
    expressionTokensDisplay(expression_tokens);
end


% Parse

% Find pairs of parentheses
parentheses_stack = [];
parentheses_ranges = zeros(0, 2);
for i = 1:length(expression_tokens)
    if strcmp(expression_tokens(i).type, 'lparen')
        parentheses_stack(end+1) = i;
    elseif strcmp(expression_tokens(i).type, 'rparen')
        if length(parentheses_stack) == 0
            error('expression cannot be fully evaluated: unmatched '')'' (right parenthesis)');
        end
        parentheses_ranges(end+1, :) = [parentheses_stack(end), i];
        parentheses_stack = parentheses_stack(1:end-1);
    end
end
if length(parentheses_stack) ~= 0
    error('expression cannot be fully evaluated: unmatched ''('' (left parenthesis)');
end

if options.verbose
    fprintf('    parentheses_stack =''%s''\n', mat2str(parentheses_stack));
    fprintf('    parentheses_ranges =''%s''\n', mat2str(parentheses_ranges));
end


% Process first to last
token_string_types_info = struct('name', {}, 'type', {}, 'operator_token', {}, 'operator_function', {});
token_string_types_info(end+1) = struct('name', 'power', 'type', 'binary_right', 'operator_token', '^', 'operator_function', @power);
token_string_types_info(end+1) = struct('name', 'divide', 'type', 'binary_left', 'operator_token', '/', 'operator_function', @rdivide);
token_string_types_info(end+1) = struct('name', 'times', 'type', 'binary', 'operator_token', '*', 'operator_function', @times);
token_string_types_info(end+1) = struct('name', 'minus', 'type', 'binary_left', 'operator_token', '-', 'operator_function', @minus);
token_string_types_info(end+1) = struct('name', 'plus', 'type', 'binary', 'operator_token', '+', 'operator_function', @plus);
token_string_types_info(end+1) = struct('name', 'parentheses', 'type', 'parentheses', 'operator_token', {{}}, 'operator_function', {{}});
% token_string_types_info

% warning('Associativity not fully checked');

while length(expression_tokens) > 1
    previous_n_expression_tokens = length(expression_tokens);
    
    for i = 1:length(token_string_types_info)
        token_string_type_info = token_string_types_info(i);
        token_string_name = token_string_type_info.name;
        token_string_type = token_string_type_info.type;
        token_string_operator_token = token_string_type_info.operator_token;
        token_string_operator_function = token_string_type_info.operator_function;
        
        is_binary = any(strcmp(token_string_type, {'binary', 'binary_left', 'binary_right'}));
        is_parentheses = strcmp(token_string_type, 'parentheses');
        is_right_associative = strcmp(token_string_type, 'binary_right');
        
        if is_binary || is_parentheses
            if is_right_associative
                j_values = length(expression_tokens)-1:-1:2;
            else
                j_values = 2:length(expression_tokens)-1;
            end
            j_values_reduced = j_values;
            while length(j_values_reduced) > 0
                j_reduced = j_values_reduced(1);
                previous_expression_token = expression_tokens(j_reduced-1);
                expression_token = expression_tokens(j_reduced);
                next_expression_token = expression_tokens(j_reduced+1);
                previous_expression_token_type = previous_expression_token.type;
                expression_token_type = expression_token.type;
                next_expression_token_type = next_expression_token.type;
                previous_expression_token_token = previous_expression_token.token;
                expression_token_token = expression_token.token;
                next_expression_token_token = next_expression_token.token;
                
                if is_binary
                    if ~strcmp(expression_token_token, token_string_operator_token) || ~any(strcmp(previous_expression_token_type, {'float', 'variable'})) || ~any(strcmp(next_expression_token_type, {'float', 'variable'}))
                        j_values_reduced = j_values_reduced(2:end);
                        continue;
                    end
                    % if ischar(previous_expression_token_token)
                        previous_expression_token_token = valueToTokenValue(previous_expression_token_token);
                    % end
                    % if ischar(next_expression_token_token)
                        next_expression_token_token = valueToTokenValue(next_expression_token_token);
                    % end
                    if options.verbose
                        fprintf('    previous_expression_token_token = %s\n', valueToString(previous_expression_token_token));
                        fprintf('    next_expression_token_token = %s\n', valueToString(next_expression_token_token));
                        fprintf('    token_string_operator_function result = %s\n', valueToString(token_string_operator_function(previous_expression_token_token, next_expression_token_token)));
                    end
                    expression_token = struct('type', 'float', 'token', token_string_operator_function(previous_expression_token_token, next_expression_token_token));
                    
                elseif is_parentheses
                    if ~strcmp(previous_expression_token_type, 'lparen') || ~strcmp(next_expression_token_type, 'rparen')
                        j_values_reduced = j_values_reduced(2:end);
                        continue;
                    end
                    if ischar(expression_token_token)
                        expression_token_token = valueToTokenValue(expression_token_token);
                    end
                    % expression_token = struct('type', 'float', 'token', expression_token_token);
                end
                
                expression_tokens = [expression_tokens(1:j_reduced-2), expression_token, expression_tokens(j_reduced+2:end)];
                j_values_reduced = j_values_reduced(3:end);
                if ~is_right_associative
                    j_values_reduced = j_values_reduced - 2;
                end
            end
            % error('Unfinished');
            if options.verbose
                expressionTokensDisplay(expression_tokens);
            end
        else
            error('Unfinished');
        end
    end
    
    if previous_n_expression_tokens == length(expression_tokens)
        % error('expression cannot be fully evaluated: no more token strings can be evaluated');
        break
    end
end


% Create final expression string:
expression = [expression_tokens(:).token];

value = valueToReturnValue(expression);


end
