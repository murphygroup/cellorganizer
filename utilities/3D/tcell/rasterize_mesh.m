function [result] = rasterize_mesh(given_mesh, options)
  % 2013-04-01 tebuck: Created.
  %
  % Test:
  % addpath(genpath('Mesh_voxelisation'))
  % test_mesh = struct;
  % test_mesh.vertices = [0, 0, 0; 1, 0, 0; 0, 1, 0; 0, 0, 1] * 20 + 1;
  % test_mesh.faces = [1, 3, 2; 1, 3, 4; 1, 4, 2; 4, 3, 2];
  % imshow(reshape_contrast(rasterize_mesh(test_mesh, struct('oversampling_scale', 4))))
  % imshow(reshape_contrast(rasterize_mesh(test_mesh, struct('oversampling_scale', 4, 'cropping', {[1, 25; 1, 30; 1, 35]}))))
  %
  % Dependencies:
  % From the File Exchange: Mesh_voxelisation
  % Murphy lab: my image_resize_nd
  %
  
  default_options.cropping = 'tight';
  default_options.oversampling_scale = 2;

  if ~exist('options', 'var')
    options = default_options;
  else
    options = process_options_structure(default_options, options);
  end

  % Is the scaled rounding necessary?
  % What is the proper way to choose voxels? Centers?
  given_mesh.vertices = given_mesh.vertices .* options.oversampling_scale;
  % spacing = 1 / options.oversampling_scale;
  % strcmpi(options.cropping, 'tight')
  if strcmpi(options.cropping, 'tight')
    % x_range = floor(min(given_mesh.vertices(:, 1))):ceil(max(given_mesh.vertices(:, 1)));
    % y_range = floor(min(given_mesh.vertices(:, 2))):ceil(max(given_mesh.vertices(:, 2)));
    % z_range = floor(min(given_mesh.vertices(:, 3))):ceil(max(given_mesh.vertices(:, 3)));
    x_range = floor(min(given_mesh.vertices(:, 1)) / options.oversampling_scale) * options.oversampling_scale:ceil(max(given_mesh.vertices(:, 1)) / options.oversampling_scale) * options.oversampling_scale;
    y_range = floor(min(given_mesh.vertices(:, 2)) / options.oversampling_scale) * options.oversampling_scale:ceil(max(given_mesh.vertices(:, 2)) / options.oversampling_scale) * options.oversampling_scale;
    z_range = floor(min(given_mesh.vertices(:, 3)) / options.oversampling_scale) * options.oversampling_scale:ceil(max(given_mesh.vertices(:, 3)) / options.oversampling_scale) * options.oversampling_scale;
  else
    % options.cropping = options.cropping .* options.oversampling_scale;
    options.cropping(:, 1) = floor(options.cropping(:, 1) ./ options.oversampling_scale) .* options.oversampling_scale^2;
    options.cropping(:, 2) = floor(options.cropping(:, 2) ./ options.oversampling_scale) .* options.oversampling_scale^2;
    x_range = options.cropping(1, 1):options.cropping(1, 2);
    y_range = options.cropping(2, 1):options.cropping(2, 2);
    z_range = options.cropping(3, 1):options.cropping(3, 2);
    % x_range = floor(options.cropping(1, 1) / options.oversampling_scale) * options.oversampling_scale:ceil(options.cropping(1, 2) / options.oversampling_scale) * options.oversampling_scale;
    % y_range = floor(options.cropping(2, 1) / options.oversampling_scale) * options.oversampling_scale:ceil(options.cropping(2, 2) / options.oversampling_scale) * options.oversampling_scale;
    % z_range = floor(options.cropping(3, 1) / options.oversampling_scale) * options.oversampling_scale:ceil(options.cropping(3, 2) / options.oversampling_scale) * options.oversampling_scale;
    % x_range = floor(options.cropping(1, 1) / options.oversampling_scale) * options.oversampling_scale^2:ceil(options.cropping(1, 2) / options.oversampling_scale) * options.oversampling_scale^2;
    % y_range = floor(options.cropping(2, 1) / options.oversampling_scale) * options.oversampling_scale^2:ceil(options.cropping(2, 2) / options.oversampling_scale) * options.oversampling_scale^2;
    % z_range = floor(options.cropping(3, 1) / options.oversampling_scale) * options.oversampling_scale^2:ceil(options.cropping(3, 2) / options.oversampling_scale) * options.oversampling_scale^2;
  end

  mesh_image = voxelise(x_range, y_range, z_range, given_mesh);
  mesh_image = permute(mesh_image, [2, 1, 3]);
  
  mesh_image = image_resize_nd(double(mesh_image), 1 / options.oversampling_scale, 'bilinear');

  result = mesh_image;

end