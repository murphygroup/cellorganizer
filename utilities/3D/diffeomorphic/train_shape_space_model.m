function [shape_space_model, shape_space_approx] = train_shape_space_model(options)
% Build a shape space model from images. Takes as input only an options structure and returns a model structure. This function can be called with the same options by multiple jobs to do certain parts in parallel.
% Assumptions:
% Images are the same size and are square in XY.
% 
% Options you must specify for good results:
% number_images, the number of images used to train this model
% image_function, a function that takes an index from 1 to options.number_images and returns the corresponding image
% convergence_absolute_error, the sum of the absolute differences between any two images at which interpolation convergence will be assumed. Specifically, during an interpolation between a source and a target image, the images will be iteratively deformed towards each other. After each iteration, this convergence test will be done.
% voxel_size
% shape_aspect_ratio, the 
% save_location
% 
% Options you should specify for completeness:
% maximum_shape_space_dimension, the maximum dimensionality for the shape space to be constructed (maximum value of 7)
% explained_variance_threshold, the minimum amount of variance that must be explained by the automatically selected dimensionality of the shape space, but this can only produce a dimensionality up to maximum_shape_space_dimension (between 0 and 1)
% 
% Returns a structure containing:
% positions, the positions of each shape in the constructed shape space
% convex_hull, the convex hull of those positions
% tessellation, the Delaunay tessellation of those positions
% 
% Other options can be seen in the default_options structure at the beginning of this function.
% 
% Dependencies from File Exchange:
% DataHash
% Dependencies from the Murphy and Rohde labs:
% 
% History:
% 2013-02-13 tebuck: Copied from test_walk_ternary025.m.
% 
% To do:
% Use MDS approximations.

% 9/22/13 grj - improved readability


    % We require certain options:
    required_options = {'image_function', 'convergence_absolute_error', 'voxel_size', 'shape_aspect_ratio', 'save_location', 'number_images'};
    for index = 1:length(required_options)
        required_option = required_options{index};
        if ~isfield(options, required_option)
          error('Option %s is required', required_option)
        end
    end

    first_image = options.image_function(1);

    [M, N, P] = size(first_image);
%     [X_0, Y_0, Z_0] = meshgrid(1:N, 1:M, 1:P);
%     identity_map = {X_0, Y_0, Z_0}; 
%     clear X_0 Y_0 Z_0
    window_size = [M, N, P];

    % Default options go in a structure passed to process_options_structure:
    default_options = struct();
    % Parameters for windowed LDDMM (these should probably be in a substructure so field names do not overlap (not that they should as of 2013-02-14)):
    default_options.filter_radius = 32;
    % By default use windows the size of the images:
    default_options.window_radius = window_size(1);
    % Flat cells:
    default_options.kernel_z_radius = default_options.filter_radius * options.shape_aspect_ratio(3);
    % Flat cells:
    default_options.maximum_deformation_per_step = [1, 1, .5];
    default_options.maximum_iterations = 1000;
    default_options.just_compute_distance = true;
    % Options specific to this function:
    
    default_options.embedding = struct();
    % Combine default and given options:
    options = ml_initparam(options, default_options);
    
      % The proportion of the variance that must be explained by the shape space is used to choose the dimensionality (up to options.maximum_shape_space_dimension):
    options.embedding = ml_initparam(options.embedding, struct( ...
                            'method', 'most_complete', ...
                            'force_positive_definiteness', true, ... %% this is only for the embed_partial_distance_matrix function
                            'explained_variance_threshold', 0.9, ...
                            'maximum_shape_space_dimension', 50 ...
                            ));
                        
                        
    shape_space_filename = [options.save_location filesep 'shape_space'];
    
    
    if ~strcmpi(options.embedding.method, 'most_complete')
        options.matCompletionFunctionString = 'matrix_completion_random_columns(distances_sampled, dist_nans, options)';
    end
    
    [distances_incomplete, distances_sampled] = compute_distance_matrix(options.image_function, options.number_images, options);
    
    
    shape_space_model = embed_distance_matrix(distances_incomplete, options);
    
    shape_space_model.imfunc = options.image_function;
    shape_space_model.shape_space_options =     rmfield(options, {'image_function'});
    
    shape_space_model.numimgs = options.number_images;
%     shape_space_model.imsize = options.imsize;
    shape_space_model.name = 'diffeomorphic model of the cell';
    shape_space_model.type = 'diffeomorphic';

    shape_space_model.version = 1.0;
end
