function [varargout] = align_synapse_in_segmentation_image(segmentation_image, synapse_point_2d, options)
% function [realigned_segmentation, transformed_images] = align_synapse_in_segmentation_image(segmentation_image, synapse_point_2d, options)
% Find the orientation of the synapse in segmentation_image using one of several methods, then rotate segmentation_image so that the synapse is facing upwards.
% synapse_point_2d is in BK's YX order (i.e., RC, not HV, and he calls it XY).
%
% Dependencies:
%
% History:
% 2012-10-31 tebuck: Copied from fit_plane_to_synapse_in_segmentation_image.m.
% 2013-02-04 tebuck: If segmentation_image is a mesh structure, treat it appropriately.
% 2013-02-28 tebuck: Making sure that all of my processing takes into account that I use XYZ coordinates (X horizontal/second dimension in Matlab), while Snake3D and affine use YXZ (X vertical/first dimension in Matlab), and the input is YX.
% 2013-03-02 tebuck: Removed special processing for segmentation_image being an image instead of a mesh and just convert it to a mesh at the beginning for simplicity.

  % default_options.curvature_smoothing = 2;
  default_options.debug_plot = false;
  default_options.images_to_transform = {};
  
  if ~exist('options', 'var')
    options = default_options;
  else
    options = process_options_structure(default_options, options);
  end
  
  % Smooth shape to remove noise, e.g., from watershed on noisy image:
  % sigma = sqrt(3) * 2.;
  % sigma = sqrt(3) * 1.;
  % sigma = sqrt(3) * .75;
  sigma = sqrt(3) * .5;
  % sigma = sqrt(3) * .25;
  % sigma = 0.;
  
  % Convert segmentation_image to a mesh if it is an image:
  if ~isstruct(segmentation_image)
    smoothed_segmentation_image = segmentation_image;
    if sigma > 0
      smoothed_segmentation_image = smooth3(smoothed_segmentation_image, 'gaussian', ceil(sigma) * 6 + 1, sigma);
    end
    segmentation_image = isosurface(smoothed_segmentation_image, .5);
  end
  
  % Prepare the mesh for alignment:
  single_cell_mesh = segmentation_image;
  % Assumes just one segmentation for now:
  synapse_point_3d = synapse_point_2d;
  
  weights = pdist2(synapse_point_3d, single_cell_mesh.vertices(:, 1:2));
  % Some of these might become inf:
  % weights = 1 ./ (weights ./ mean(weights));
  weights = exp(-weights ./ mean(weights));
  weights = weights ./ sum(weights);
  synapse_point_3d(3) = weights * single_cell_mesh.vertices(:, 3);
  % synapse_point_3d

  
  % % Rotate the mesh based on the synapse position:
  % mean_vertex = mean()

  % Angle the synapse should be facing (0 is to the right, angle increases counterclockwise in images):
  desired_angle_degrees = -90;
  
  % Fit a plane to the mesh:
  % disable_initial_alignment = false;
  disable_initial_alignment = true;
  if disable_initial_alignment
    fprintf('>>>> HACK, disable_initial_alignment = true\n')
    initial_transform = eye(4);
    initial_transform_affine = eye(4);
  else
    % Rotate the mesh based on the synapse position:
    mean_vertex = mean(single_cell_mesh.vertices, 1);
    direction = synapse_point_3d - mean_vertex;
    direction = direction ./ norm(direction);
    initial_r = sqrt(direction(1)^2 + direction(2)^2 + direction(3)^2);
    initial_theta = acos(direction(3) / initial_r);
    initial_phi = atan2(direction(2), direction(1));
    initial_theta = initial_theta - pi / 2;
    initial_theta_degrees = initial_theta * 180 / pi;
    initial_phi_degrees = initial_phi * 180 / pi;
    
    % fprintf('>>>> HACK, initial_theta_degrees = 0\n'), initial_theta_degrees = 0;
    % fprintf('>>>> HACK, initial_phi_degrees = 0\n'), initial_phi_degrees = 0;
    
    initial_translation_component = synapse_point_3d;
    initial_transform = translation_matrix3h(initial_translation_component) * rotation_matrix3h(3, -initial_phi_degrees + desired_angle_degrees) * rotation_matrix3h(2, -initial_theta_degrees) * translation_matrix3h(-(initial_translation_component));
    initial_transform_affine = translation_matrix3h(initial_translation_component - 1) * rotation_matrix3h(3, -initial_phi_degrees + desired_angle_degrees) * rotation_matrix3h(2, -initial_theta_degrees) * translation_matrix3h(-(initial_translation_component - 1));
    % fprintf('>>>> HACK, initial_transform''s rotation_scale_component = eye(4)\n'), initial_transform = translation_matrix3h(initial_translation_component) * eye(4) * translation_matrix3h(-(initial_translation_component)); initial_transform_affine = translation_matrix3h(initial_translation_component - 1) * eye(4) * translation_matrix3h(-(initial_translation_component - 1));
    single_cell_mesh.vertices = [single_cell_mesh.vertices, ones(size(single_cell_mesh.vertices, 1), 1)] * initial_transform';
    single_cell_mesh.vertices = single_cell_mesh.vertices(:, 1:3);
  end

  % % Have an option to select the method for this:
  % switch options.synapse_detection_method
    % case 'curvature local minima'
      % [plane_parameters] = fit_plane_to_synapse_in_mesh(single_cell_mesh, synapse_point_3d, options);
      [synapse_parameters] = fit_plane_to_synapse_in_mesh(single_cell_mesh, synapse_point_3d, options);
      % plane_parameters = synapse_parameters.plane_parameters;
    % otherwise
      % error('Unknown synapse_detection_method value ''%s''', synapse_detection_method)
  % end

  
  % % Around which point should we rotate? Use the template's parameters? synapse_point_2d and the Z centroid of the image?
  default_template_options = get_default_template_options();
  % center = [default_template_options.yc, default_template_options.xc, default_template_options.zc];

  % Rotate to align the plane's normal:
  
  % Using the template image size from BK:
  % template_center = [default_template_options.xc, default_template_options.yc, default_template_options.zc];
  % Use XYZ order for my coordinates here:
  template_center = [default_template_options.yc, default_template_options.xc, default_template_options.zc];
  % Using one of the images' size if it exists, but should have an option for this?:
  if numel(options.images_to_transform) > 0
    template_image_size = size(options.images_to_transform{1});
    % template_center = (template_image_size - 1) * .5 + 1;
  end
  % template_image_size = [default_template_options.imx, default_template_options.imy, default_template_options.imz];
  template_array_size = [default_template_options.imy, default_template_options.imx, default_template_options.imz];

  % synapse_parameters
  n_x = synapse_parameters.synapse_normal(1);
  n_y = synapse_parameters.synapse_normal(2);
  n_z = synapse_parameters.synapse_normal(3);
  translation_component = synapse_parameters.synapse_center;
  % scale_component = (default_template_options.yc + default_template_options.zc) / 2 / synapse_parameters.synapse_radius;
  scale_component = (default_template_options.yr + default_template_options.zr) / 2 / synapse_parameters.synapse_radius;
  % default_template_options, template_image_size, pause
  % Using the template image size from BK:
  centering_component = template_center - synapse_parameters.synapse_center;
  % % Just center it in this image size:
  % synapse_parameters.synapse_center
  % size(segmentation_image)
  % keyboard
  % centering_component = size(segmentation_image) / 2 - synapse_parameters.synapse_center;

  % synapse_parameters
  % scale_component
  % r = sqrt(n_x^2 + n_y^2 + n_z^2);
  % n_z / r
  % theta = acos(n_z / r)
  % phi = atan2(n_y, n_x)
  % theta = theta - pi / 2;
  [phi, theta, r] = cart2sph(n_x, n_y, n_z);
  theta_degrees = theta * 180 / pi;
  phi_degrees = phi * 180 / pi;
  % scale_component
  % fprintf('>>>> HACK, scale_component = 1\n'), scale_component = 1;
  % fprintf('>>>> HACK, theta_degrees = 0\n'), theta_degrees = 0;
  % fprintf('>>>> HACK, phi_degrees = 0\n'), phi_degrees = 0;

  % whos smoothed_segmentation_image, pause
  % fprintf('>>>> HACK, translation_component = synapse_point_3d\n'), translation_component = synapse_point_3d;
  % fprintf('>>>> HACK, translation_component = 0\n'), translation_component = synapse_point_3d * 0;
  rotation_scale_component = rotation_matrix3h(3, -phi_degrees + desired_angle_degrees) * rotation_matrix3h(2, -theta_degrees) * scale_matrix3h(scale_component * ones(1, 3));
  % fprintf('>>>> HACK, rotation_scale_component = eye(4)\n'), rotation_scale_component = eye(4);
  transform = translation_matrix3h(centering_component) * translation_matrix3h(translation_component) * rotation_scale_component * translation_matrix3h(-(translation_component));
  transform_affine = translation_matrix3h(centering_component - 1) * translation_matrix3h(translation_component - 1) * rotation_scale_component * translation_matrix3h(-(translation_component - 1));
  % initial_transform
  
  transform = transform * initial_transform;
  transform_affine = transform_affine * initial_transform_affine;
  % Use YXZ order for affine:
  transform_affine = transform_affine([2, 1, 3, 4], [2, 1, 3, 4]);
  realigned_segmentation.vertices = [segmentation_image.vertices, ones(size(segmentation_image.vertices, 1), 1)] * transform';
  % realigned_segmentation.vertices = [segmentation_image.vertices, ones(size(segmentation_image.vertices, 1), 1)] * transform;
  realigned_segmentation.vertices = realigned_segmentation.vertices(:, 1:3);

  
  % realigned_segmentation.vertex_transform = transform;
  % realigned_segmentation.affine_input_transform = transform_affine;
  
  varargout = cell(1, nargout);
  
  if nargout >= 1
    varargout{1} = realigned_segmentation;
  end
  if nargout >= 2
    % transform
    transformed_images = cellfun(@(x)affine(x, transform_affine, [], false, [], [], 'crop'), options.images_to_transform, 'UniformOutput', false);
    % % Crop to BK template image size:
    % % template_array_size
    % % transformed_image_sizes = cell2mat(cellfun(@(x)size(x), transformed_images, 'UniformOutput', false)')
    % % transformed_images = cellfun(@(x)x(1:end, 1:end, 1:end), transformed_images, 'UniformOutput', false);
    % transformed_images = cellfun(@(x)cat(3, x(1:template_array_size(1), 1:template_array_size(2), 1:min(end, template_array_size(3))), zeros(template_array_size(1:3) - [0, 0, size(x, 3)])), transformed_images, 'UniformOutput', false);
    varargout{2} = transformed_images;
  end
  % Save the transformations for meshes and images (for the latter, assume the images will be transformed using affine from FEX, which modifies translations by 1 because of Matlab indexing):
  if nargout >= 3
    alignment_transforms = struct('vertex_transform', transform, 'affine_input_transform', transform_affine);
    varargout{3} = alignment_transforms;
  end
  
  