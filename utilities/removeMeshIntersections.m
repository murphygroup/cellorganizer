function [ meshData ] = removeMeshIntersections( meshData, min_clearance, method )
%REMOVEMESHINTERSECTIONS Adjust a set of triangle meshes to enforce clearance.
%
% Inputs
% ------
% meshData      = 
% min_clearance = 
% method        = 
%
% Outputs
% -------
% meshData = meshData with converted CSGdata appended


% Authors: Taraz Buck
%
% Copyright (C) 2020 Murphy Lab
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

if nargin < 2 || isempty(min_clearance)
    min_clearance = 0;
end

if min_clearance ~= 0
    error('CellOrganizer:instance2MCellMDL:notImplemented', 'min_clearance ~= 0 not implemented');
end

if nargin < 3 || isempty(method)
    % method = 'delete_higher_index';
    method = 'delete_intersects_framework';
end


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
meshData_flat_essential_mask = meshData_flat_boundaries_mask;
switch method
    case {'delete_higher_index'}
        meshData_flat_to_start = true(1, length(meshData_flat));
    case {'delete_intersects_framework'}
        meshData_flat_to_start = meshData_flat_boundaries_mask;
end
% TODO: For framework meshes, flatten intersecting patches and move back along normals by min_clearance
% TODO: For other meshes, use a closest point algorithm or approximate with kd-tree to enforce min_clearance

use_octree = false;

% Find self-intersections
% for i = find(~meshData_flat_domain_mask)
for i = 1:length(meshData_flat)
    item = meshData_flat(i);
    item_mesh = item.mesh;
    % Do not keep intersecting meshes
    if use_octree
        error('CellOrganizer:instance2MCellMDL:notImplemented', 'Not implemented');
        % mesh_self_intersects = fastMesh2Mesh(item_mesh.vertices, item_mesh.vertices, item_mesh.faces, item_mesh.faces, octs, true);
    else
        mesh_self_intersects = any(SurfaceIntersection(item_mesh, 'self'), 1:2);
    end
    % meshData_flat_keep_mask(i) = meshData_flat_keep_mask(i) && ~mesh_self_intersects;
    if mesh_self_intersects
        if meshData_flat_essential_mask(i)
            error('CellOrganizer:removeMeshIntersections:selfIntersection', 'Framework or domain mesh has self-intersection');
        end
        meshData_flat_keep_mask(i) = false;
    end
end

switch method
    case {'delete_higher_index', 'delete_intersects_framework'}
        % For every intersecting pair of meshes, remove the mesh with the higher index (usually a vesicle)
        
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
                    error('CellOrganizer:instance2MCellMDL:notImplemented', 'Not implemented');
                    % meshes_intersect = fastMesh2Mesh(mesh1.vertices, mesh2.vertices, mesh1.faces, mesh2.faces, octs, true);
                else
                    meshes_intersect = any(SurfaceIntersection(mesh1, mesh2), 1:2);
                end
                if meshes_intersect && meshData_flat_essential_mask(j)
                    error('CellOrganizer:removeMeshIntersections:intersection', 'Framework or domain meshes intersect');
                end
                meshData_flat_keep_mask(j) = ~meshes_intersect;
            end
        end
        
        k = length(meshData_flat_keep_mask);
        for i = length(meshData):-1:1
            for j = length(meshData(i).list):-1:1
                if ~meshData_flat_keep_mask(k)
                    meshData(i).list = [meshData(i).list(1:j - 1), meshData(i).list(j + 1:end)];
                end
                k = k - 1;
            end
        end
        
        n_removed = sum(~meshData_flat_keep_mask);
        fprintf('Removed %d intersecting meshes\n', n_removed);
        
    otherwise
        error('CellOrganizer:removeMeshIntersections:unknownMethod', 'Unknown method %s', method);
end


end
