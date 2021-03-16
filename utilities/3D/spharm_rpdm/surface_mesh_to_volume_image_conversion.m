function [bim] = surface_mesh_to_volume_image_conversion(input_mesh, img_size, options)
% convert surface mesh to volume image (binary)
% The function will not check whether the input mesh is within the range of
% image size, so make sure the surfaces are places within the proposed
% image. 

% xruan 05/09/2019 add check of mesh size versus img_size, if img_size is
% small, then throw out a warning

addpath(genpath('Mesh_voxelisation'));

if nargin < 2
    img_size = max(vertices);
end

default_options.cropping = 'tight';
default_options.oversampling_scale = 2;

if ~exist('options', 'var')
    options = default_options;
else
    options = process_options_structure(default_options, options);
end

rasterization_oversampling_scale = options.oversampling_scale;
segmentation_rasterization_function = @(given_mesh, window_size)rasterize_mesh(given_mesh, struct('oversampling_scale', rasterization_oversampling_scale, 'cropping', [ones(3, 1), window_size']));

% input_mesh.vertices = Zvert_pdm;
% input_mesh.faces = fs;

% compar mesh size and img_size
mesh_max_size = max(input_mesh.vertices);
if any(mesh_max_size > img_size([2, 1, 3]))
    warning('The provided image size [%0.0f %0.0f %0.0f] is smaller than the mesh size [%0.1f %0.1f %0.1f] for at least one dimension! You may consider use large image size by change options.model.spharm_rpdm.imageSize!', ...
        img_size(2), img_size(1), img_size(3), mesh_max_size(1), mesh_max_size(2), mesh_max_size(3));
end

bim = segmentation_rasterization_function(input_mesh, img_size([2, 1, 3]));
% bim = imresize3(int8(bim > 0.5), [256, 256, 60], 'nearest');

if any(size(bim) > img_size)
    bim = crop_image_function(bim, img_size);
elseif any(size(bim) < img_size)
    bim = pad_image_3d_function(bim, img_size);    
end
    

if false
    figure, patch('faces', input_mesh.faces, 'vertices', input_mesh.vertices, 'FaceVertexCData',jet(size(input_mesh.vertices,1)), 'FaceColor', 'interp')
    view([45, 45]);
    figure, imshow(reshape_2d(bim > 0.5, -1), [])
end


end