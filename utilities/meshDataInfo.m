function [ info ] = meshDataInfo( meshData )
%MESHDATAINFO Adjust a set of triangle meshes to enforce clearance.
%
% Inputs
% ------
% meshData = 
%
% Outputs
% -------
% info = meshData flattened, masks, etc.


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

nucleus_mesh_index = find(meshData_flat_framework_mask & strcmp(meshData_flat_subnames, 'NU'));
cell_mesh_index = find(meshData_flat_framework_mask & strcmp(meshData_flat_subnames, 'CP'));
ec_mesh_index = find(meshData_flat_framework_mask & strcmp(meshData_flat_subnames, 'EC'));

info = struct();
info.meshData_flat = meshData_flat;
info.meshData_flat_names = meshData_flat_names;
info.meshData_flat_subnames = meshData_flat_subnames;
info.meshData_flat_framework_mask = meshData_flat_framework_mask;
info.meshData_flat_domain_mask = meshData_flat_domain_mask;
info.meshData_flat_boundaries_mask = meshData_flat_boundaries_mask;
info.nucleus_mesh_index = nucleus_mesh_index;
info.cell_mesh_index = cell_mesh_index;
info.ec_mesh_index = ec_mesh_index;


end
