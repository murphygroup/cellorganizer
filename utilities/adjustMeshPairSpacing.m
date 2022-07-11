function [ mesh1_mod, mesh2_mod, meshes_intersect, meshes_n_intersections ] = adjustMeshPairSpacing( mesh1, mesh2, options_output )
%ADJUSTMESHPAIRSPACING Adjust a set of triangle meshes to enforce clearance.
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
% Copyright (C) 2022 Murphy Lab
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

min_clearance = options_output.framework_min_clearance;
intersecting_mesh_object_policy = options_output.intersecting_mesh_object_policy;
maximum_filter_rounds_value = options_output.framework_clearance_n_max_filter_rounds;

separation_distribution = 1;
nucleus_normal_proportion = 1;
maximum_filter_rounds_value = 1;
if isfield(options_output, 'separation_distribution')
    separation_distribution = options_output.separation_distribution;
end
if isfield(options_output, 'nucleus_normal_proportion')
    nucleus_normal_proportion = options_output.nucleus_normal_proportion;
end
if isfield(options_output, 'maximum_filter_rounds_value')
    maximum_filter_rounds_value = options_output.maximum_filter_rounds_value;
end

if min_clearance > -inf
    % Move cell vertices along nucleus normals in an attempt to remove intersections (e.g., nucleus sticking out of cell)
    
    hypotn = @(x)sqrt(sum(x.^2, 2));
    scalar_product = @(x, y)sum(x .* y, 2);
    normalize = @(x)x ./ hypotn(x);
    
    % tr = triangulation(mesh1.faces, mesh1.vertices)
    mesh1_adjacency = logical(meshAdjacencyMatrix(mesh1.faces));
    % mesh2_adjacency = logical(meshAdjacencyMatrix(mesh2.faces));
    
    n1 = vertexNormal(mesh1.vertices, mesh1.faces);
    v1 = mesh1.vertices;
    v2 = mesh2.vertices;
    v1_mean = mean(v1, 1);
    v1 = v1 - v1_mean;
    v2 = v2 - v1_mean;
    v2_minus_v1 = v2 - v1;
    
    v2_minus_v1_normalized = normalize(v2_minus_v1);
    v1_normalized = normalize(v1);
    v2_minus_v1_angles = acos(max(-1, min(1, sum(v2_minus_v1_normalized .* v1_normalized, 2))));
    
    v2_minus_v1_projection_v1 = scalar_product(v2_minus_v1, v1_normalized);
    v2_minus_v1_projection_n1 = scalar_product(v2_minus_v1, n1);
    
    
    % Impose minimum projections
    modification_projection_v1 = -min(0, v2_minus_v1_projection_v1 - min_clearance);
    modification_projection_n1 = -min(0, v2_minus_v1_projection_n1 - min_clearance);
    % Spread values by maximum filter to reduce self-intersection that would otherwise remain
    for i = 1:maximum_filter_rounds_value
        temp = modification_projection_v1;
        for j = 1:size(mesh1.vertices, 1)
            temp(j) = max(temp(mesh1_adjacency(j, :)));
        end
        modification_projection_v1 = temp;
        temp = modification_projection_n1;
        for j = 1:size(mesh1.vertices, 1)
            temp(j) = max(temp(mesh1_adjacency(j, :)));
        end
        modification_projection_n1 = temp;
    end
    % Add modified projection to cell and nucleus meshes
    modification = (1 - nucleus_normal_proportion) * modification_projection_v1 .* v1_normalized;
    modification = modification + nucleus_normal_proportion * modification_projection_n1 .* n1;
    v1_modification = -(1 - separation_distribution) * modification;
    v2_modification = separation_distribution * modification;
    v1_mod = v1 + v1_modification;
    v2_mod = v2 + v2_modification;
    
    v1_mod_normalized = normalize(v1_mod);
    v2_mod_minus_v1_mod = v2_mod - v1_mod;
    v2_mod_minus_v1_mod_projection_v1_mod = scalar_product(v2_mod_minus_v1_mod, v1_mod_normalized);
    
    mesh1_mod = mesh1;
    mesh1_mod.vertices = v1_mod + v1_mean;
    mesh2_mod = mesh2;
    mesh2_mod.vertices = v2_mod + v1_mean;
    
    meshes_intersect = [];
    meshes_n_intersections = [];
    if nargout > 2
        meshes_n_intersections = sum(SurfaceIntersection(mesh1_mod, mesh2_mod), 1:2);
        % meshes_intersect = any(SurfaceIntersection(mesh1_mod, mesh2_mod), 1:2);
        meshes_intersect = meshes_n_intersections > 0;
    end
else
    mesh1_mod = mesh1;
    mesh2_mod = mesh2;
end


end
