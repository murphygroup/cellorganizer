function [vertices_out, faces_out] = reorder_mesh_vertices(vertices, faces, method, direction)
% the function is used for reorder mesh vertices, so that in the patch, the
% color scheme will be the same

if nargin < 3
    method = 'z_axis';
end

nface = size(faces, 1);
nvert = size(vertices, 1);
orig_inds = 1 : nvert;

switch method
    case 'x-axis'
        [vertices_out, idx] = sortrows(vertices, [1, 2], {'descend', 'descend'});
    case 'y-axis'
        [vertices_out, idx] = sortrows(vertices, [2, 3], {'descend', 'descend'});
    case 'z-axis'
        [vertices_out, idx] = sortrows(vertices, [3, 1], {'descend', 'descend'});
    case 'given-direction'
        direction = direction / norm(direction);
        direction_1 = direction([2, 1, 3]) .* [1, -1, 0];
        if all(direction_1 == 0)
            % in this case, the direction is z-axis
            direction_1 = [1, 0, 0]; 
        end
        direction_mat = [direction; direction_1];
        project_vertices = vertices * direction_mat';
        [~, idx] = sortrows(project_vertices, [1, 2], {'descend', 'descend'});     
        vertices_out = vertices(idx, :);
end

[~, idx_1] = sort(idx);       
mapper = containers.Map(orig_inds, orig_inds(idx_1));
faces_out = arrayfun(@(x) mapper(x), faces(:));
faces_out = reshape(faces_out, nface, []);

  
end
