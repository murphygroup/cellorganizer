function [ result ] = double2str( value )
%DOUBLE2STR Convert double to string losslessly.

if isnumeric(value)
    % result = sprintf('%f', value);
    value = sprintf('%.17e', value);
end

if ischar(value)
    value = strtrim(value);
    if any(strcmpi(value, {'nan', 'inf', '+inf', '-inf'}))
        result = value;
    else
        value2 = str2double(strtrim(value));
        
        if ~isnan(value2)
            result = value;
        else
            error('Argument must be numeric or char representing float');
        end
    end
else
    error('Argument must be numeric or char representing float');
end

% Remove trailing zeros
result_tokens = regexp(result, '(\.)([0-9]*?)(0*)(e)', 'tokens');
if length(result_tokens) > 0
    result_tokens = result_tokens{1};
    old_string = strjoin(result_tokens, '');
    if isempty(result_tokens{2})
        result_tokens{1} = '';
    end
    result_tokens{3} = '';
    new_string = strjoin(result_tokens, '');
    result = strrep(result, old_string, new_string);
end

result_tokens = regexp(result, '(e[+-])(0*)([0-9]*?)$', 'tokens');
if length(result_tokens) > 0
    result_tokens = result_tokens{1};
    old_string = strjoin(result_tokens, '');
    if isempty(result_tokens{3})
        result_tokens = {};
    else
        result_tokens{2} = '';
    end
    new_string = strjoin(result_tokens, '');
    result = strrep(result, old_string, new_string);
end

end