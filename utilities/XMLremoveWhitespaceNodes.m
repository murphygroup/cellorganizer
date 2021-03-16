function [node] = XMLremoveWhitespaceNodes(node)
% XMLREMOVEWHITESPACENODES Remove whitespace text nodes

% Copyright (C) 2018 Taraz Buck
% Copyright (C) 2019 Murphy Lab
% Computational Biology Department
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

% 2019-02-18 Taraz Buck: Copied `setupVCML.m` to `XMLremoveWhitespaceNodes.m`.

node_list = {node};
child_index_list = [0];
nodes_to_remove = {};
nodes_to_remove_parents = {};
deleted_a_child = false;
while length(node_list) > 0
    if ~deleted_a_child
        child_list = node_list{end}.getChildNodes();
        if child_list.getLength() > 0
            % Go down a level
            node_list{end+1} = child_list.item(child_index_list(end));
            child_index_list(end+1) = 0;
            continue;
        end
    end
    deleted_a_child = false;
    if strcmp(char(class(node_list{end})), 'org.apache.xerces.dom.DeferredTextImpl') ...
            && length(strtrim(char(node_list{end}.getTextContent()))) == 0
        % Remove blank text node
        nodes_to_remove{end+1} = node_list{end};
        nodes_to_remove_parents{end+1} = node_list{end-1};
    end
    child_index_list(end) = child_index_list(end) + 1;
    if child_index_list(end) >= node_list{end}.getChildNodes().getLength()
        % Go up a level
        node_list = node_list(1:end-1);
        child_index_list = child_index_list(1:end-1);
        % Go to next child
        deleted_a_child = true;
        continue;
    end
end
for i = 1:length(nodes_to_remove)
    node_to_remove = nodes_to_remove{i};
    % node_to_remove.getParentNode().removeChild(node_to_remove);
    node_to_remove_parent = nodes_to_remove_parents{i};
    node_to_remove_parent.removeChild(node_to_remove);
end

