function value = mathEval(expression, variables_map, options)
%MATHEVAL Evaluates arithmetic expression string without EVAL.
% VALUE = MATHEVAL(EXPRESSION) returns the double value of the evaluated string EXPRESSION or throws an error indicating failure.
%
% Inputs
% ------
% expression    = string containing an arithmetic expression
% variables_map = optional containers.Map from identifier strings to expression strings
% options       = optional struct containing the following optional fields (default values in parentheses):
%     identifier_pattern ('[a-zA-Z_][a-zA-Z0-9_]*')
%     expression_max_length (1e4)
%     return_double (true)
%     verbose (false)
% 
% Outputs
% -------
% value = double value of expression


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


if nargin < 2
    variables_map = containers.Map('KeyType', 'char', 'ValueType', 'char');
end

if nargin < 3
    options = struct();
end

% Process options
default_options = struct();
% Variable name
default_options.identifier_pattern = '[a-zA-Z_][a-zA-Z0-9_]*';
% Maximum length of transformed expression before throwing an error (protect against recursive variable expressions)
default_options.expression_max_length = 1e4;
default_options.return_double = true;
default_options.verbose = false;

if ~exist('options', 'var')
    options = default_options; 
else
    options = process_options_structure(default_options, options);
end

% Postprocess options
identifier_pattern = ['(?:',options.identifier_pattern,')'];


ws = '(?:[ \t]+)';
wsopt = '(?:[ \t]*)';

operators_patterns = {};
operators_patterns{end+1} = '+';
operators_patterns{end+1} = '\-';
operators_patterns{end+1} = '*';
operators_patterns{end+1} = '/';
operators_patterns{end+1} = '\^';

operators_pattern = ['[',cell2mat(reshape(operators_patterns, 1, [])),']'];

function x = isolatePattern(x)
    x = ['(?<=',operators_pattern,wsopt,'|\(',wsopt,'|^)',x,'(?=',wsopt,operators_pattern,'|',wsopt,'\)|$)'];
end

% int_pattern = '(?:[0-9]+)';
% TODO: Edit this to match an authoritative grammar
float_pattern = ['(?:[\-\+]?(?:(?:(?:(?<![0-9.])[0-9][0-9]*(?![0-9.])|[0-9]+\.[0-9]*(?![0-9])|(?<![0-9])[0-9]*\.[0-9]+)(?:[eE][\-\+]?[0-9]+)?)|[iI][nN][fF]|[nN][aA][nN]))'];
% Test: regexp('5 5. .5 0.5 15 15. .15 15.15 -5 -5. -.5 -0.5 -15 -15. -.15 -15.15 5e3 5.e3 .5e3 0.5e3 15e3 15.e3 .15e3 15.15e3 +5e-3 +5.e-3 +.5e-3 +0.5e-3 +15e-3 +15.e-3 +.15e-3 +15.15e-3 inf Inf INF nan NaN NAN', '(?:[\-\+]?(?:(?:(?:(?<![0-9.])[0-9][0-9]*(?![0-9.])|[0-9]+\.[0-9]*(?![0-9])|(?<![0-9])[0-9]*\.[0-9]+)(?:[eE][\-\+]?[0-9]+)?)|[iI][nN][fF]|[nN][aA][nN]))', 'match')
float_pattern = isolatePattern(float_pattern);

% identifier_or_float_pattern = ['(?:',identifier_pattern,'|',float_pattern,')'];

% variable_pattern = ['(?<variable_name>(?<=[+\-*/]',wsopt,'|\(',wsopt,'|^)',identifier_pattern,'(?=',wsopt,'[+\-*/]|',wsopt,'\)|$))'];
% variable_pattern = ['(?<variable_name>(?<=',operators_pattern,wsopt,'|\(',wsopt,'|^)',identifier_pattern,'(?=',wsopt,operators_pattern,'|',wsopt,'\)|$))'];
variable_pattern = ['(?<variable_name>',isolatePattern(identifier_pattern),')'];

% Process first to last
token_types_info = struct('type', {}, 'pattern', {});
token_types_info(end+1) = struct('type', 'lparen', 'pattern', '\(');
token_types_info(end+1) = struct('type', 'rparen', 'pattern', '\)');
token_types_info(end+1) = struct('type', 'float', 'pattern', float_pattern);
token_types_info(end+1) = struct('type', 'power', 'pattern', '\^');
token_types_info(end+1) = struct('type', 'divide', 'pattern', '/');
token_types_info(end+1) = struct('type', 'times', 'pattern', '\*');
token_types_info(end+1) = struct('type', 'minus', 'pattern', '\-');
token_types_info(end+1) = struct('type', 'plus', 'pattern', '\+');
token_types_info(end+1) = struct('type', 'variable', 'pattern', variable_pattern);


function x = isStringNumeric(x)
    x = ~isnan(str2double(x));
end

function x = num2str_double(x)
    x = num2str(x, 18);
end


% Tokenize

expression_tokens = struct('type', {}, 'token', {});
expression_unprocessed_portion = expression;
expression_length_processed = 0;
while length(expression_unprocessed_portion) > 0
    token_found = false;
    for i = 1:length(token_types_info)
        token_type_info = token_types_info(i);
        token_type = token_type_info.type;
        token_pattern = ['^', token_type_info.pattern];
        tokens = regexp(expression_unprocessed_portion, token_pattern, 'match');
        if length(tokens) == 0
            continue;
        end
        token = tokens{1};
        expression_tokens(end+1) = struct('type', token_type, 'token', token);
        expression_length_processed = expression_length_processed + length(token);
        if expression_length_processed > options.expression_max_length
            error('expression cannot be fully evaluated: transformed expression too long');
        end
        expression_unprocessed_portion = expression_unprocessed_portion(length(token)+1:end);
        
        if strcmp(token_type, 'variable')
            if ~variables_map.isKey(token)
                error('expression cannot be fully evaluated: unknown variable ''%s''', token);
            end
            expression_unprocessed_portion = [variables_map(token), expression_unprocessed_portion];
            expression_tokens = expression_tokens(1:end-1);
        end
        token_found = true;
        break;
    end
    if ~token_found
        error('expression cannot be fully evaluated: unrecognized token beginning with ''%s''', expression_unprocessed_portion(1:10));
    end
end

if options.verbose
    expression_tokens_cell = [fieldnames(expression_tokens)'; squeeze(struct2cell(expression_tokens))'];
    fprintf('    expression =\n'); disp(expression)
    fprintf('    expression_tokens_cell =\n'); disp(expression_tokens_cell)
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
                    if ~strcmp(expression_token_token, token_string_operator_token) || ~strcmp(previous_expression_token_type, 'float') || ~strcmp(next_expression_token_type, 'float')
                        j_values_reduced = j_values_reduced(2:end);
                        continue;
                    end
                    if ischar(previous_expression_token_token)
                        previous_expression_token_token = str2double(previous_expression_token_token);
                    end
                    if ischar(next_expression_token_token)
                        next_expression_token_token = str2double(next_expression_token_token);
                    end
                    expression_token = struct('type', 'float', 'token', token_string_operator_function(previous_expression_token_token, next_expression_token_token));
                    
                elseif is_parentheses
                    if ~strcmp(previous_expression_token_type, 'lparen') || ~strcmp(next_expression_token_type, 'rparen')
                        j_values_reduced = j_values_reduced(2:end);
                        continue;
                    end
                    if ischar(expression_token_token)
                        expression_token_token = str2double(expression_token_token);
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
                % ['type', 'token'; {expression_tokens.type}', {expression_tokens.token}']
                expression_tokens_cell = [fieldnames(expression_tokens)'; squeeze(struct2cell(expression_tokens))'];
                fprintf('    expression_tokens_cell =\n'); disp(expression_tokens_cell)
            end
        else
            error('Unfinished');
        end
    end
    
    if previous_n_expression_tokens == length(expression_tokens)
        error('expression cannot be fully evaluated: no more token strings can be evaluated');
    end
end


if length(expression_tokens) ~= 1
    error('expression cannot be fully evaluated: final token string has length other than one');
end
if ~strcmp(expression_tokens.type, 'float')
    error('expression cannot be fully evaluated: final token is not of type ''float''');
end

expression = expression_tokens.token;
is_expression_numeric_string = ischar(expression) && isStringNumeric(expression);
if ~isfloat(expression) && ~is_expression_numeric_string
    error('expression cannot be fully evaluated: final transformed expression not numeric');
end

value = expression;
if options.return_double && is_expression_numeric_string
    value = str2double(value);
elseif ~options.return_double
    value = num2str(value);
end


end
