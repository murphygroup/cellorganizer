function [ CSGdata, meshData ] = convertCSGToMesh( CSGdata, meshData, models, imgs, options )
%CONVERTCSGTOMESH Convertes CellOrganizer geometry in CSGdata format to meshData format.
%
% Inputs
% ------
% CSGdata  = 
% meshData = 
%
% Outputs
% -------
% meshData = meshData with converted CSGdata appended


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


debug = false;
% debug = true;

warning_id = ['CellOrganizer:', mfilename];
if ~debug
    warning('off', warning_id);
end

if ~isfield(options.output, 'use_individual_meshes')
    options.output.use_individual_meshes = false;
    % options.output.use_individual_meshes = true;
end
use_individual_meshes = options.output.use_individual_meshes;

% name_map = options.output.name_map;
% resolution_before_downsampling = options.output.resolution_before_downsampling;
model_names = options.output.model_names;


warning(warning_id, 'This should be done further upstream in createSBMLFrameworkstruct so meshes are scaled like objects (see createSBMLstruct3) and then adjusted to output_length_unit here');
for i = 1:length(meshData)
    for j = 1:length(meshData(i).list)
        single_meshData_list_item_resolution = meshData(i).list(j).resolution;
        % single_meshData_list_item_resolution = resolution_before_downsampling;
        if isfield(meshData(i).list(j), 'mesh')
            meshData(i).list(j).mesh.vertices = meshData(i).list(j).mesh.vertices .* repmat(single_meshData_list_item_resolution, size(meshData(i).list(j).mesh.vertices, 1), 1);
        end
        meshData(i).list(j).resolution = [1, 1, 1];
    end
end

meshData(:).source = 'mesh';
%{
for i = 1:length(meshData)
    for j = 1:length(meshData(i).list)
        meshData(i).list(j).name = translateWithDefaultIdentity(name_map, meshData(i).list(j).name);
    end
end
%}

CSGdata_fieldnames = fieldnames(CSGdata);
empty_mesh = struct('vertices', zeros(0, 3), 'faces', zeros(0, 3));
CSGdata_list_lengths = zeros(length(CSGdata_fieldnames));
for i = 1:length(CSGdata_fieldnames)
    CSGdata_fieldname = CSGdata_fieldnames{i};
    single_CSGdata = CSGdata.(CSGdata_fieldname);
    CSGdata_list_lengths(i) = length(single_CSGdata.list);
end
CSGdata_list_lengths_cumsum = [0, cumsum(CSGdata_list_lengths(1:end-1))];
    
for i = 1:length(CSGdata_fieldnames)
    CSGdata_fieldname = CSGdata_fieldnames{i};
    single_CSGdata = CSGdata.(CSGdata_fieldname);
    object_name = '';
    object_name_is_EC = strcmp(CSGdata_fieldname, 'EC');
    
    meshData_name = 'modelMesh';
    if object_name_is_EC
        meshData_name = 'domainMesh';
    end
    
    single_meshData = struct('name', meshData_name, 'source', 'CSG', 'list', struct('type', {}, 'name', {}, 'ordinal', {}, 'mesh', {}, 'resolution', {}, 'img', {}));
    
    object_img = [];
    switch CSGdata_fieldname
        case 'EC'
            object_name = CSGdata_fieldname;
        case 'vesicle'
            if length(models) > 1
                warning(warning_id, 'More than one protein model given. Assuming all vesicles are from first protein model.');
            end
            object_name = model_names{1};
            switch lower(options.synthesis)
                case {'all', 'framework'}
                    imgs_offset = 2;
                case {'cell', 'nucleus'}
                    imgs_offset = 1;
            end
            object_img = imgs{imgs_offset+1};
        otherwise
            warning(warning_id, 'CSGdata contains unrecognized field "%s"', CSGdata_fieldname);
    end
    
    if isempty(object_name)
        continue;
    end
    if ~use_individual_meshes
        single_CSGdata_mesh = empty_mesh;
    end
    for j = 1:length(single_CSGdata.list)
        single_CSGdata_list_item = single_CSGdata.list(j);
        single_CSGdata_list_item_resolution = single_CSGdata_list_item.resolution;
        item_mesh = empty_mesh;
        if ~isfield(single_CSGdata_list_item, 'scale')
            single_CSGdata_list_item.position = ones(1, 3);
        end
        if ~isfield(single_CSGdata_list_item, 'rotationmatrix')
            single_CSGdata_list_item.rotationmatrix = eye(4);
        end
        if ~isfield(single_CSGdata_list_item, 'position')
            single_CSGdata_list_item.position = [0, 0, 0];
        end
        switch single_CSGdata_list_item.type
            % case {'cube', 'sphere'}
            case 'cube'
                item_mesh = createCube();
                warning(warning_id, 'Check cube');
                item_mesh.vertices = item_mesh.vertices - 0.5;
                item_mesh.vertices = item_mesh.vertices * 2;
                item_mesh = rmfield(item_mesh, 'edges');
            % case 'sphere (disabled)'
            case 'sphere'
                %{
                item_mesh = struct();
                [item_mesh.vertices, item_mesh.faces] = createIcosahedron();
                item_mesh = SubdivideSphericalMesh(item_mesh, 2);
                %}
                %{
                item_mesh = createTetrahedron();
                item_mesh.vertices = item_mesh.vertices - 0.5;
                item_mesh.vertices = item_mesh.vertices * sqrt(4/3);
                item_mesh = SubdivideSphericalMesh(item_mesh, 3);
                %}
                item_mesh = struct();
                [item_mesh.vertices, item_mesh.faces] = SpiralSampleSphere(10 * 4^2 + 2);
            otherwise
                warning(warning_id, 'CSGdata contains unrecognized field "%s"', CSGdata_fieldname);
        end
        
        if ~isempty(item_mesh.vertices)
            if size(item_mesh.faces,2) == 4
                item_mesh.faces = [item_mesh.faces(:, 1:3); item_mesh.faces(:, [1, 3, 4])];
            end
            % fprintf('*** convertCSGToMesh: item_mesh.vertices mean = %s, std = %s\n', num2str(mean(item_mesh.vertices, 1)), num2str(std(item_mesh.vertices, 0, 1)));
            item_mesh.vertices = item_mesh.vertices .* repmat(single_CSGdata_list_item.scale, size(item_mesh.vertices, 1), 1);
            item_mesh.vertices = item_mesh.vertices * single_CSGdata_list_item.rotationmatrix(1:3, 1:3)';
            item_mesh.vertices = item_mesh.vertices + repmat(single_CSGdata_list_item.position, size(item_mesh.vertices, 1), 1);
            % warning(warning_id, 'item_mesh.vertices not translated"');
            
            item_mesh.vertices = item_mesh.vertices .* repmat(single_CSGdata_list_item_resolution, size(item_mesh.vertices, 1), 1);
            if options.output.SBMLFlipXToAlign
                % Flip in X to align with image output
                max_size_voxels = size(options.cell);
                max_size_voxels_xyz = max_size_voxels([2, 1, 3]);
                FV = item_mesh;
                FV.vertices(:, 1) = (max_size_voxels_xyz(1) - (FV.vertices(:, 1) ./ options.resolution.cell(1) - 1)) .* options.resolution.cell(1);
                FV.faces(:, 2:3) = FV.faces(:, [3, 2]);
                item_mesh = FV;
            end
            
            if use_individual_meshes
                single_meshData_list_item = struct('type', 'triangle mesh', 'name', object_name, 'ordinal', CSGdata_list_lengths_cumsum(i) + j, 'mesh', item_mesh, 'resolution', [1, 1, 1], 'img', []);
                single_meshData.list(end+1) = single_meshData_list_item;
            else
                item_mesh.faces = item_mesh.faces + size(single_CSGdata_mesh.vertices, 1);
                single_CSGdata_mesh.vertices = [single_CSGdata_mesh.vertices; item_mesh.vertices];
                single_CSGdata_mesh.faces = [single_CSGdata_mesh.faces; item_mesh.faces];
            end
        end
    end
    
    if use_individual_meshes && ~isempty(single_meshData.list)
        single_meshData.list(1).img = object_img;
        meshData(end+1) = orderfields(single_meshData, meshData);
    elseif ~use_individual_meshes && ~isempty(single_CSGdata_mesh.vertices)
        %{
        object_name = translateWithDefaultIdentity(name_map, object_name);
        % if ~is_membrane_function(object_name)
        if network_info.compartments(object_name).spatial_dimensions == 2
            object_name = network_info.compartments(object_name).outside;
        end
        %}
        compartment_object_name = object_name;
        compartment_mesh_name = compartment_object_name;
        % compartment_mesh_name = [compartment_object_name, '_mesh'];
        
        single_meshData_list_item = struct('type', 'triangle mesh', 'name', compartment_mesh_name, 'ordinal', i, 'mesh', single_CSGdata_mesh, 'resolution', [1, 1, 1], 'img', object_img);
        single_meshData.list(end+1) = single_meshData_list_item;
        meshData(end+1) = orderfields(single_meshData, meshData);
    end
end


end
