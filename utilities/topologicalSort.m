function [ keys_topological_order ] = topologicalSort( objects, parents_function )
%TOPOLOGICALSORT Sort nodes in a directed acyclic graph.
%
% Inputs
% ------
% objects          = array, cell array, or struct array of objects
% parents_function = handle to function that returns the parents of an object
%
% Outputs
% -------
% keys_topological_order = array of indices (if objects is an array or cell array) or cell array of keys (if objects is a container.Map)
%
% Example
% -------
% objects = [3, 1, 2]; parents_function = @(x)num2cell(objects(objects < x)); topologicalSort(objects, parents_function)


% Authors: Taraz Buck
%
% Copyright (C) 2019 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
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


% add_parentless_keys_immediately = false;
add_parentless_keys_immediately = true;

% Produce topological sort of objects (Kahn 1962, "Topological sorting of large networks")
if isnumeric(objects) || iscell(objects)
    keys = num2cell(1:length(objects))';
elseif isa(objects, class(containers.Map))
    keys = objects.keys()';
else
    error('objects is not of class array, cell array, or struct array');
end
keys_topological_order = {};
if length(keys) == 0
    return;
end
% Initialize 
parentless_keys = {};
keys_children = containers.Map('KeyType', class(keys{1}), 'ValueType', 'any');
keys_parents = containers.Map('KeyType', class(keys{1}), 'ValueType', 'any');
for i = 1:length(keys)
    key = keys{i};
    key_parents = parents_function(key);
    if length(key_parents) == 0
        parentless_keys{end+1} = key;
    else
        keys_parents(key) = key_parents;
        for j = 1:length(key_parents)
            keys_parent = key_parents{j};
            if ~keys_children.isKey(keys_parent)
                keys_children(keys_parent) = {};
            end
            keys_children(keys_parent) = [keys_children(keys_parent), {key}];
        end
    end
end
% Iteratively find parentless keys and add them to the end of the list
% L = keys_topological_order
% S = parentless_keys
% n = current_parentless_key
if add_parentless_keys_immediately
    keys_topological_order = [keys_topological_order, parentless_keys];
end
while ~isempty(parentless_keys)
    % Remove a parentless node from the graph and add it to the sorted list
    current_parentless_key = parentless_keys{end};
    parentless_keys(end) = [];
    if ~add_parentless_keys_immediately
        keys_topological_order{end+1} = current_parentless_key;
    end
    if ~keys_children.isKey(current_parentless_key)
        continue;
    end
    % Add any children that now have no parents to parentless_keys
    current_children_keys = keys_children(current_parentless_key);
    keys_children.remove(current_parentless_key);
    for i = 1:length(current_children_keys)
        current_child_key = current_children_keys{i};
        current_child_parents = keys_parents(current_child_key);
        for j = 1:length(current_child_parents)
            current_child_parent = current_child_parents{j};
            if objectsEqual(current_child_parent, current_parentless_key)
                current_child_parents = [current_child_parents(1:j-1), current_child_parents(j+1:end)];
                keys_parents(current_child_key) = current_child_parents;
                break;
            end
        end
        if length(current_child_parents) == 0
            keys_parents.remove(current_child_key);
            parentless_keys{end+1} = current_child_key;
            if add_parentless_keys_immediately
                keys_topological_order{end+1} = current_child_key;
            end
        end
    end
end
if ~isempty(keys_parents) || ~isempty(keys_children)
    error('Directed graph is not acyclic');
end
if length(keys_topological_order) ~= length(keys)
    error('Topological sort failed');
end


end




function value = objectsEqual(object1, object2)

value = false;
value = value || (isnumeric(object1) && isnumeric(object2) && object1 == object2);
value = value || (ischar(object1) && ischar(object2) && strcmp(object1, object2));

end

