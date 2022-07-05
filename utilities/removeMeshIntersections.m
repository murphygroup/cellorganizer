function [ meshData ] = removeMeshIntersections( meshData, options )
%REMOVEMESHINTERSECTIONS Adjust a set of triangle meshes to enforce clearance.
%
% Inputs
% ------
% meshData = 
% options  = 
%
% Outputs
% -------
% meshData = meshData with converted CSGdata appended


% Authors: Taraz Buck
%
% Copyright (C) 2020-2022 Murphy Lab
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


warning_id = ['CellOrganizer:', mfilename];
warning('off', warning_id);
warning('on', warning_id); % Debug

min_clearance = options.output.framework_min_clearance;
intersecting_mesh_object_policy = options.output.intersecting_mesh_object_policy;
maximum_filter_rounds_value = options.output.framework_clearance_n_max_filter_rounds;

check_intersections = ~strcmp(intersecting_mesh_object_policy, 'ignore');
remove_on_any_intersection = strcmp(intersecting_mesh_object_policy, 'remove');
remove_on_framework_intersection = strcmp(intersecting_mesh_object_policy, 'remove_framework');
reject_on_any_intersection = strcmp(intersecting_mesh_object_policy, 'reject');
reject_on_framework_intersection = strcmp(intersecting_mesh_object_policy, 'reject_framework')

% Put all meshes into a flat struct array for ease of processing
meshData_flat = repmat(setfield(meshData(1).list(1), 'super_name', ''), 0, 1);
meshData_flat_names = {};
meshData_flat_subnames = {};
for i = 1:length(meshData)
    for j = 1:length(meshData(i).list)
        meshData_flat_names{end + 1} = meshData(i).name;
        meshData_flat_subnames{end + 1} = meshData(i).list(j).name;
        meshData_flat(end + 1) = setfield(meshData(i).list(j), 'super_name', meshData_flat_names{end});
    end
end

meshData_flat_framework_mask = strcmp(meshData_flat_names, 'frameworkMesh');
meshData_flat_domain_mask = strcmp(meshData_flat_names, 'domainMesh');
meshData_flat_boundaries_mask = meshData_flat_framework_mask | meshData_flat_domain_mask;
meshData_flat_keep_mask = true(1, length(meshData_flat));
meshData_flat_essential_mask = false(1, length(meshData_flat));
if reject_on_framework_intersection
    meshData_flat_essential_mask = meshData_flat_essential_mask | meshData_flat_boundaries_mask;
end
if reject_on_any_intersection
    meshData_flat_essential_mask = meshData_flat_essential_mask | true;
end
meshData_flat_modified_mask = false(1, length(meshData_flat));;

meshData_flat_to_start = false(1, length(meshData_flat));
if remove_on_framework_intersection
    % delete_intersects_framework
    meshData_flat_to_start = meshData_flat_to_start | meshData_flat_boundaries_mask;
end
if remove_on_any_intersection
    % delete_higher_index
    meshData_flat_to_start = meshData_flat_to_start | true;
end

nucleus_mesh_index = find(strcmp(meshData_flat_names, 'frameworkMesh') & strcmp(meshData_flat_subnames, 'NU'));
nucleus_mesh = meshData_flat(nucleus_mesh_index).mesh;
nucleus_mesh_adjacency = logical(meshAdjacencyMatrix(nucleus_mesh.faces));
cell_mesh_index = find(strcmp(meshData_flat_names, 'frameworkMesh') & strcmp(meshData_flat_subnames, 'CP'));
cell_mesh = meshData_flat(cell_mesh_index).mesh;

options_output = options.output;
% options_output.separation_distribution = 1;
% options_output.nucleus_normal_proportion = 1;
% options_output.maximum_filter_rounds_value = 1;
[nucleus_mesh_mod, cell_mesh_mod] = adjustMeshPairSpacing(nucleus_mesh, cell_mesh, options_output);

meshData_flat(nucleus_mesh_index).mesh = nucleus_mesh_mod;
meshData_flat(cell_mesh_index).mesh = cell_mesh_mod;

%{
if min_clearance > -inf
    % Move cell vertices along nucleus normals in an attempt to remove intersections (nucleus sticking out of cell)
    
    nn = vertexNormal(nucleus_mesh.vertices, nucleus_mesh.faces);
    nv = nucleus_mesh.vertices;
    cv = cell_mesh.vertices;
    nv_mean = mean(nv, 1);
    nv = nv - nv_mean;
    cv = cv - nv_mean;
    cv_minus_nv = cv - nv;
    
    cv_minus_nv_projection_nn = scalar_product(cv_minus_nv, nn);

    % Impose minimum projections
    modification_projection_nn = -min(0, cv_minus_nv_projection_nn - min_clearance);
    % Spread values by maximum filter to reduce self-intersection that would otherwise remain
    for i = 1:maximum_filter_rounds_value
        temp = modification_projection_nn;
        for j = 1:size(nucleus_mesh.vertices, 1)
            temp(j) = max(temp(nucleus_mesh_adjacency(j, :)));
        end
        modification_projection_nn = temp;
    end
    % Add modified projection to cell and nucleus meshes
    modification = modification_projection_nn .* nn;
    cv_modification = separation_distribution * modification;
    cv_mod = cv + cv_modification;
    
    cell_mesh.vertices = cv_mod + nv_mean;
    
    meshData{cell_mesh_index} = adjustMeshPairSpacing(nucleus_mesh, cell_mesh, options);
end
%}


% TODO: For other meshes, use a closest point algorithm or approximate with kd-tree to enforce min_clearance

use_octree = false;

if check_intersections
    % Find self-intersections
    for i = 1:length(meshData_flat)
        item = meshData_flat(i);
        item_mesh = item.mesh;
        % Do not keep intersecting meshes
        if use_octree
            error('CellOrganizer:removeMeshIntersections:notImplemented', 'Not implemented');
            % mesh_self_intersects = fastMesh2Mesh(item_mesh.vertices, item_mesh.vertices, item_mesh.faces, item_mesh.faces, octs, true);
        else
            % mesh_self_intersects = any(SurfaceIntersection(item_mesh, 'self'), 1:2);
            mesh_self_intersects = any(SurfaceIntersectionBlock(item_mesh, 'self', 4), 1:2);
        end
        if mesh_self_intersects
            if meshData_flat_essential_mask(i)
                error('CellOrganizer:removeMeshIntersections:selfIntersection', 'Framework or domain mesh has self-intersection');
            end
            meshData_flat_keep_mask(i) = false;
        end
    end
end

% For every intersecting pair of meshes, remove the mesh with the higher index (usually a vesicle)

if check_intersections
    % Maximum bin size for fastMesh2Mesh's octree
    octs = 10;
    % Compare all pairs of meshes
    for i = find(meshData_flat_to_start)
        name1 = meshData_flat(i).name;
        % Keep framework and EC
        if ~meshData_flat_keep_mask(i)
            continue
        end
        mesh1 = meshData_flat(i).mesh;
        for j = i + 1:length(meshData_flat)
            name2 = meshData_flat(j).name;
            if meshData_flat_essential_mask(j) || ~meshData_flat_keep_mask(j)
                continue
            end
            mesh2 = meshData_flat(j).mesh;
            
            % Do not keep intersecting meshes
            if use_octree
                error('CellOrganizer:removeMeshIntersectionsremoveMeshIntersections:notImplemented', 'Not implemented');
                % meshes_intersect = fastMesh2Mesh(mesh1.vertices, mesh2.vertices, mesh1.faces, mesh2.faces, octs, true);
            else
                % meshes_intersect = any(SurfaceIntersection(mesh1, mesh2), 1:2);
                meshes_intersect = any(SurfaceIntersectionBlock(mesh1, mesh2, 4), 1:2);
            end
            if meshes_intersect && meshData_flat_essential_mask(j)
                error('CellOrganizer:removeMeshIntersections:intersection', 'Framework or domain meshes intersect');
            end
            meshData_flat_keep_mask(j) = ~meshes_intersect;
        end
    end
end

k = length(meshData_flat_keep_mask);
for i = length(meshData):-1:1
    for j = length(meshData(i).list):-1:1
        if ~meshData_flat_keep_mask(k)
            meshData(i).list = [meshData(i).list(1:j - 1), meshData(i).list(j + 1:end)];
        elseif meshData_flat_modified_mask(k)
            meshData(i).list(j).mesh = meshData_flat(k).mesh;
        end
        k = k - 1;
    end
end

n_removed = sum(~meshData_flat_keep_mask);
n_modified = sum(~meshData_flat_keep_mask & meshData_flat_modified_mask);
fprintf('Removed %d intersecting meshes, modified %d\n', n_removed, n_modified);


end
