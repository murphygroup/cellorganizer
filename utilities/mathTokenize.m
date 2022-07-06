function expression_tokens = mathTokenize(expression, variables_map, options)
%MATHTOKENIZE Tokenizes arithmetic expression string without EVAL.
% EXPRESSION_TOKENS = MATHTOKENIZE(EXPRESSION) returns a structure representing the tokenization of the string EXPRESSION or throws an error indicating failure.
%
% Inputs
% ------
% expression    = char containing an arithmetic expression
% options       = optional struct containing the following optional fields (default values in parentheses):
%     identifier_pattern ('[a-zA-Z_][a-zA-Z0-9_]*')
%     expression_max_length (1e4)
%     check_variables (true)
%     verbose (false)
% 
% Outputs
% -------
% expression_tokens = struct array containing the following optional fields (types in parentheses):
%     type (char)
%     token (char)


% Copyright (C) 2019 Taraz Buck
% Copyright (C) 2019 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
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


% 2019-07-16 Taraz Buck: Created.


if nargin < 2 || isempty(variables_map)
    variables_map = containers.Map('KeyType', 'char', 'ValueType', 'char');
end

if nargin < 3
% if nargin < 2
    options = struct();
end

% Process options
default_options = struct();
% Variable name
default_options.identifier_pattern = '[a-zA-Z_][a-zA-Z0-9_]*';
% Maximum length of transformed expression before throwing an error (protect against recursive variable expressions)
default_options.expression_max_length = 1e4;
default_options.check_variables = true;
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
% float_pattern = ['(?:[\-\+]?(?:(?:(?:(?<![0-9.])[0-9][0-9]*(?![0-9.])|[0-9]+\.[0-9]*(?![0-9])|(?<![0-9])[0-9]*\.[0-9]+)(?:[eE][\-\+]?[0-9]+)?)|[iI][nN][fF]|[nN][aA][nN]))'];
float_pattern = '(?:(?<![0-9.])(?:[\-\+]?(?:(?:(?:[0-9]+|[0-9]+\.[0-9]*|[0-9]*\.[0-9]+)(?:[eE][\-\+]?[0-9]+)?)|[iI][nN][fF]|[nN][aA][nN]))(?![0-9.]))';
% Test: regexp('5 5. .5 0.5 15 15. .15 15.15 -5 -5. -.5 -0.5 -15 -15. -.15 -15.15 5e3 5.e3 .5e3 0.5e3 15e3 15.e3 .15e3 15.15e3 +5e-3 +5.e-3 +.5e-3 +0.5e-3 +15e-3 +15.e-3 +.15e-3 +15.15e-3 inf Inf INF nan NaN NAN', float_pattern, 'match')
float_pattern = isolatePattern(float_pattern);

% identifier_or_float_pattern = ['(?:',identifier_pattern,'|',float_pattern,')'];

% variable_pattern = ['(?<variable_name>(?<=[+\-*/]',wsopt,'|\(',wsopt,'|^)',identifier_pattern,'(?=',wsopt,'[+\-*/]|',wsopt,'\)|$))'];
% variable_pattern = ['(?<variable_name>(?<=',operators_pattern,wsopt,'|\(',wsopt,'|^)',identifier_pattern,'(?=',wsopt,operators_pattern,'|',wsopt,'\)|$))'];
variable_pattern = ['(?<variable_name>',isolatePattern(identifier_pattern),')'];

% Process first to last
token_types_info = struct('type', {}, 'pattern', {});
token_types_info(end+1) = struct('type', 'ignore', 'pattern', ws);
token_types_info(end+1) = struct('type', 'lparen', 'pattern', '\(');
token_types_info(end+1) = struct('type', 'rparen', 'pattern', '\)');
% token_types_info(end+1) = struct('type', 'float', 'pattern', float_pattern);
token_types_info(end+1) = struct('type', 'power', 'pattern', '\^');
token_types_info(end+1) = struct('type', 'divide', 'pattern', '/');
token_types_info(end+1) = struct('type', 'times', 'pattern', '\*');
token_types_info(end+1) = struct('type', 'minus', 'pattern', '\-');
token_types_info(end+1) = struct('type', 'plus', 'pattern', '\+');
% token_types_info(end+1) = struct('type', 'ignore', 'pattern', ws);
token_types_info(end+1) = struct('type', 'float', 'pattern', float_pattern);
token_types_info(end+1) = struct('type', 'variable', 'pattern', variable_pattern);



% Enable use of multiple variable types including numeric, string, and DimensionedExpression variables

function x = num2str_double(x)
    x = num2str(x, 18);
end

variables_map_keys = variables_map.keys();
if length(variables_map_keys) > 0
    variable_type = class(variables_map(variables_map_keys{1}));
else
    variable_type = 'char';
end
switch variable_type
    case 'double'
        valueToString = @(x)double2str(x);
        stringToToken = @(x)str2double(x);
        concatenateTokens = @(varargin)cat(varargin{:});
    case 'DimensionedExpression'
        valueToString = @(x)x.value;
        %{
        % valueToString = @(x)x;
        valueToString = @(x)char(x);
        %}
        stringToToken = @(x)DimensionedExpression(x);
        concatenateTokens = @(varargin)prod(varargin{:});
    case 'char'
        valueToString = @(x)x;
        stringToToken = @(x)x;
        % concatenateTokens = @(varargin)cat(varargin{:});
        concatenateTokens = @(varargin)strcat(varargin{:});
end

%{
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
%}

%{
function x = concatenateTokens(varargin)
    if length(varargin) == 0
        x = '';
    else
        switch class(varargin{1})
            case 'double'
                x = double2str(x);
            case 'DimensionedExpression'
                x = char(x);
            case 'char'
                x = x;
        end
    end
end
%}


% Tokenize

expression_tokens = struct('type', {}, 'token', {});
expression_unprocessed_portion = expression;
expression_length_processed = 0;
if options.verbose
    fprintf('    expression = %s\n', expression);
end
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
        
        if strcmp(token_type, 'variable') && options.check_variables
            if ~variables_map.isKey(token)
                error('expression cannot be fully evaluated: unknown variable ''%s''', token);
            end
            % Replace and continue
            if options.verbose
                fprintf('    replacing ''%s'' with ''%s''\n', token, valueToString(variables_map(token)));
            end
            expression_unprocessed_portion = ['(', valueToString(variables_map(token)), ')',  expression_unprocessed_portion];
            expression_tokens = expression_tokens(1:end-1);
        end
        
        token_found = true;
        break;
    end
    if ~token_found
        error('expression cannot be fully evaluated: unrecognized token beginning with ''%s''', expression_unprocessed_portion(1:min(10, end)));
    end
    if strcmp(token_type, 'ignore')
        expression_tokens = expression_tokens(1:end-1);
    end
    if options.verbose
        expression_tokens_cell = [fieldnames(expression_tokens)'; squeeze(struct2cell(expression_tokens))'];
        fprintf('    expression_tokens_cell =\n'); disp(expression_tokens_cell)
    end
end

%{
% Replace variables

function x = isStringNumeric(x)
    x = ~isnan(str2double(x));
end

expression_tokens2 = struct('type', {}, 'token', {});
for i = 1:length(expression_tokens)
    expression_token = expression_tokens(i);
    token_type = expression_token.type;
    token = expression_token.token;
    if ~strcmp(token_type, 'variable')
        while variables_map.isKey(expression_token)
            token = variables_map(token);
        end
        if ischar(token)
            if isStringNumeric(token)
                expression_token.type = 'float';
                expression_token.token = str2double(token);
            else
                expression_token = mathTokenize(expression, variables_map, options);
            end
        elseif isa(token, 'DimensionedExpression')
            % if token.hasVariables()
            if ~isStringNumeric(expression.value)
                expression_token = mathTokenize(expression.value, variables_map, options);
            end
        end
    end
    expression_tokens2(end+1:end+length(expression_token)) = expression_token;
end
%}

% Apply signs
expression_tokens2 = repmat(expression_tokens(1), 1, 0);
for i = 1:length(expression_tokens)
    if i > 2 && i <= length(expression_tokens) && any(strcmp(expression_tokens(i-2).type, {'lparen', 'power', 'divide', 'times', 'minus', 'plus'})) && any(strcmp(expression_tokens(i-1).type, {'minus', 'plus'})) && strcmp(expression_tokens(i).type, 'float')
        expression_tokens2(end) = expression_tokens(i);
        % expression_tokens2(end).token = [expression_tokens(end-1).token, expression_tokens(end).token];
        % expression_tokens2(end).token = concatenateTokens(expression_tokens(end-1).token, expression_tokens(end).token);
        switch expression_tokens(i-1).type
            case 'plus'
                expression_tokens2(end).token = expression_tokens(end).token;
            case 'minus'
                expression_tokens2(end).token = -expression_tokens(end).token;
        end
    else
        expression_tokens2(end+1) = expression_tokens(i);
    end
end
expression_tokens = expression_tokens2;

if options.verbose
    expression_tokens_cell = [fieldnames(expression_tokens)'; squeeze(struct2cell(expression_tokens))'];
    fprintf('    expression = %s\n', expression);
    fprintf('    expression_tokens_cell =\n'); disp(expression_tokens_cell)
end


end
