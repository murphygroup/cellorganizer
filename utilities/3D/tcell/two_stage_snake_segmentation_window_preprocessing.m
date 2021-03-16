function [preprocessing_structure] = two_stage_snake_segmentation_window_preprocessing(given_image, options)
  % Preprocess for two_stage_snake_segmentation. Not very generic.
  % 
  % 
  % 2014-04-26 tebuck: Copied from two_stage_snake_segmentation.m.
  % 
  % 
  % 
  % 
  % 
  % % Dependencies:
  % % From the File Exchange: export_fig, pmkmp, toolbox_graph, toolbox_fast_marching, affine, Snake3D, Mesh_voxelisation, dirr, inhull
  % % From tebuck: modified version of Jieyue Li's HPA_lib, rasterize_mesh

  
  default_options = struct();
  % default_options.background_subtraction_function = @(x)x;
  default_options.background_subtraction_function = @(x)x - min(x(:));
  default_options.window_center_2d = [];
  default_options.cell_diameter = [];
  
  
  if ~exist('options', 'var')
    options = default_options;
  else
    options_prior_to_processing = options;
    options = process_options_structure(default_options, options);
  end
  
  
  
  % Compute a few modified version of the image and other metadata useful for visualization and segmentation:
  image_to_segment = given_image;

  % Background subtraction (prevents edges at image boundaries):
  image_to_segment = options.background_subtraction_function(image_to_segment);
  % Scaling for consistency in force computations:
  image_to_segment = image_to_segment ./ prctile(image_to_segment(:), 99.9);
  % Hopefully makes gradients at higher intensities less influential:
  image_to_segment = reshape(histeq(reshape_2d(image_to_segment), 64), size(image_to_segment));

  preprocessing_structure = struct();
  preprocessing_structure.raw_image = given_image;
  % preprocessing_structure.background_subtracted_image = background_subtracted_image;
  preprocessing_structure.image_to_segment = image_to_segment;
  % Anisotropic filtering:
  preprocessing_structure.coherent_image_to_segment = CoherenceFilter(preprocessing_structure.image_to_segment, struct('T', 50, 'dt', 2, 'Scheme', 'R', 'verbose', 'none'));
  [gx, gy, gz] = tcell_smooth_gradient(preprocessing_structure.coherent_image_to_segment, 'scharr5', false, true);
  preprocessing_structure.fine_white_line_image_to_segment = sqrt(gx.^2 + gy.^2 + gz.^2);
  preprocessing_structure.coarse_white_line_image_to_segment = smooth3(preprocessing_structure.fine_white_line_image_to_segment, 'gaussian', 7, sqrt(3)).^2;
  
  preprocessing_structure.synapse_locations_3d = [options.window_center_2d, size(given_image, 3) / 2 + 1];
  centroid = image_centroid(reshape(image_to_segment(round(preprocessing_structure.synapse_locations_3d(2)), round(preprocessing_structure.synapse_locations_3d(1)), :), 1, []));
  preprocessing_structure.synapse_locations_3d(3) = centroid(1);
  sigma = options.cell_diameter * .2;
  smoothed_image_to_segment = convnfft_fast(preprocessing_structure.image_to_segment, gaussian_kernel(sigma));
  preprocessing_structure.refined_synapse_locations_3d = climb_image(smoothed_image_to_segment, preprocessing_structure.synapse_locations_3d);
    
  if false
    % Debug plot and such:
    % deconvolved = deconvblind(preprocessing_structure.raw_image, ones(size(preprocessing_structure.raw_image)));
    % deconvolved = deconvblind(preprocessing_structure.raw_image, ones([1, 1, 1] * options.cell_diameter));
    % deconvolved = deconvblind(preprocessing_structure.raw_image, ones([2, 2, 1] * options.cell_diameter));
    % close
    % imshow(reshape_contrast([preprocessing_structure.raw_image; deconvolved], -1));
    % a = CoherenceFilter(preprocessing_structure.raw_image, struct('T', 50, 'dt', 2, 'Scheme', 'R', 'verbose', 'none'));
    % a = CoherenceFilter(preprocessing_structure.raw_image, struct('T', 50, 'dt', 2, 'Scheme', 'R', 'verbose', 'none', 'sigma', 2));
    a = CoherenceFilter(preprocessing_structure.raw_image, struct('T', 50, 'dt', 2, 'Scheme', 'R', 'verbose', 'none', 'sigma', .25));
    [deconvolved, psf] = deconvblind(a, ones(size(preprocessing_structure.raw_image)));
    close
    imshow(reshape_contrast([a; deconvolved], -1));
    imshow(reshape_contrast(psf, -1));
    keyboard
  end
  
  
end

