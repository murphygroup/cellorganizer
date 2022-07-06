function value = translateWithDefaultIdentity(map_object, key, suffixes, recursion_limit)
%TRANSLATEWITHDEFAULTIDENTITY Recursively translates KEY using MAP_OBJECT
%
% Inputs
% ------
% map_object = 
% key        = 
%
% Outputs
% -------
% value = 


% Authors: Taraz Buck
%
% Copyright (C) 2019-2020 Murphy Lab


if nargin < 3
    suffixes = {};
end
if nargin < 4
    recursion_limit = 20;
end

% Preprocess suffixes
if ischar(suffixes)
    if size(suffixes, 1) > 1
        suffixes = arrayfun(@(x)suffixes(x, :), 1:size(suffixes, 1), 'UniformOutput', false);
    else
        suffixes = {suffixes};
    end
end
suffixes = sort(suffixes);
suffixes_preprocessed = reverse(sortrows(reverse(reshape(suffixes, [], 1))));

value = key;

% Detect suffix
suffixes_preprocessed_matches_indexes = strends(value, suffixes_preprocessed);
suffixes_preprocessed_match_index = find(suffixes_preprocessed_matches_indexes, 1, 'last');
if ~isempty(suffixes_preprocessed_match_index)
    suffixes_preprocessed_match = suffixes_preprocessed{suffixes_preprocessed_match_index};
    value = value(1:end - length(suffixes_preprocessed_match));
end

previous_value = '';
i = 0;
while ~strcmp(value, previous_value) && map_object.isKey(value) && i <= recursion_limit
    value = map_object(value);
    i = i + 1;
end

if ~isempty(suffixes_preprocessed_match_index)
    value = [value, suffixes_preprocessed_match];
end

end
