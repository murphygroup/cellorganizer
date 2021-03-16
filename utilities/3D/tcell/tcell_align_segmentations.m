function [t_cell_info] = tcell_align_segmentations(t_cell_info, options)
% The inputs are t_cell_info, options, the former one contains the general 
% information for the pipeline, and the latter contains the information for 
% the specific cell. The output is  t_cell_info. And the output result relative 
% to the cell is saved in the disk.
% 
% The idea is to use the segmented image in the disk, and align the cell so 
% that the synapse region approximately face up. We use the left and right 
% end point of the synapse region from the annotation and the centroid of the 
% cell as landmarks. The z coordinates of the left and right end points are
% inferred by weighted intensity in the neighborhood. The we used procrustes 
% analysis with the corresponding landmarks in the template to get the 
% transformation matrix.  Further more, we rescaled the image to the volume 
% of the template. And the image is then transformed based on the transformation 
% matrix. And the aligned image is saved as a binary image in the disk.
%
% 2016-02-29 xruan: Copied from master_script_align_segmentations_new.m.
% 2016-02-21 xruan: change the parallel computing. If the aligned mat
% exists, then just skip it but not save as an empty file. 
% 
% 
% 
% % Dependencies:
% % From the File Exchange: export_fig, pmkmp, toolbox_graph, toolbox_fast_marching, affine, Snake3D, Mesh_voxelisation, dirr, inhull
% % From tebuck: modified version of Jieyue Li's HPA_lib, rasterize_mesh
%
% Author: Xiongtao Ruan

  
  % Variables referenced more than once can be copied to local variables:
  master_script_options = t_cell_info.options;
  regions_location = t_cell_info.path_info.regions_location;
  segmentations_filtered_location = t_cell_info.path_info.segmentations_filtered_location;
  alignments_location = t_cell_info.path_info.alignments_location;
  deformations_location = t_cell_info.path_info.deformations_location;
  
  template_centroid = t_cell_info.template_info.template_centroid;
  template_synapse = t_cell_info.template_info.template_synapse;
  template_synapse_segment = t_cell_info.template_info.template_synapse_segment;
  template_image = t_cell_info.template_info.template_image;

  window_size_2d = t_cell_info.preprocessing_info.window_size_2d;
  window_center_2d = t_cell_info.preprocessing_info.window_center_2d;
  cropped_size = t_cell_info.preprocessing_info.cropped_size;
    
  segmentation_rasterization_function = t_cell_info.segmentation_info.segmentation_rasterization_function;
  
  template_volume = sum(template_image(:));
  % This should eventually be larger to eliminate degenerate segmentations, e.g., 5% of template volume and/or high surface-to-volume ratio:
  % segmentation_volume_threshold = 0;
  % Must be more than options.segmentation_volume_minimum_ratio of template_volume:
  segmentation_volume_threshold = template_volume * master_script_options.segmentation_volume_minimum_ratio;
  
  image_name = options.image_name;
  synapse_annotation = options.synapse_annotation;
  relative_time = options.relative_time;
  frame_channel = options.frame_channel;
  frame_index = options.frame_index;
  run_index = frame_index(:, 2);  

  % Could use either of these possibly faster methods than used in segment_using_snakes:
  % 1. Only align frames pair-wise.
  % 2. Only align landmarks.

  debug_print_alignment_skip_reasons = true;
  
  function [result_landmark] = segmentation_2d_to_3d_landmark_function(given_segmentation_image, given_landmark)
    % Returns a segmentation-derived 3D annotation from a 2D annotation by a weighted centroid in Z.
    [x, y, z] = meshgrid(1:size(given_segmentation_image, 2), 1:size(given_segmentation_image, 1), 1:size(given_segmentation_image, 3));
    horizontal_distances = sqrt((x - given_landmark(1)).^2 + (y - given_landmark(2)).^2);
    % Normalize distances:
    median_distance = median(horizontal_distances(given_segmentation_image(:) > 0));
    distance_scale = 1;
    horizontal_weights = exp(-(horizontal_distances .* distance_scale).^2 ./ median_distance);
    % if false
      % % Debug plot:
      % imshow(contrast_stretch(mean(horizontal_weights, 3))), pause
      % imshow(contrast_stretch(mean(given_segmentation_image, 3))), pause
      % imshow(contrast_stretch(mean(horizontal_weights .* given_segmentation_image, 3))), pause
    % end
    % segmentation_volume = sum(given_segmentation_image(:));
    horizontal_weights = horizontal_weights .* given_segmentation_image;
    horizontal_weights = horizontal_weights ./ sum(horizontal_weights(:));
    synapse_z = sum(z(:) .* horizontal_weights(:));
    
    centroid = image_centroid(given_segmentation_image);
    
    result_landmark = [given_landmark, synapse_z];
  end
  
  function [result_landmark] = segmentation_2d_to_3d_landmark_function_mean(given_segmentation_image, given_landmark)
    % Returns a segmentation-derived 3D annotation from a 2D annotation by a weighted centroid in Z.
    slice_num = size(given_segmentation_image, 3);
    synapse_z = sum(reshape(sum(sum(given_segmentation_image, 1), 2), [1, slice_num]) .* [1 : slice_num]) ...
        / sum(given_segmentation_image(:), 1);
    
    centroid = image_centroid(given_segmentation_image);
    
    result_landmark = [given_landmark, synapse_z];
  end


  function [result_landmark] = segmentation_2d_to_3d_landmark_function_max(given_segmentation_image, given_landmark)
    % Returns a segmentation-derived 3D annotation from a 2D annotation by a weighted centroid in Z.
    total_intensity = sum(sum(given_segmentation_image, 1), 2);
    [~, synapse_z] = max(total_intensity); 
    result_landmark = [given_landmark, synapse_z];
  end


  function [result_landmarks] = segmentation_landmark_function(given_window_image, given_segmentation_image, given_annotations, align_method)
    % Returns synapse and segmentation-derived centroid.
    if nargin < 4
        align_method = 'original'
    end
    if master_script_options.use_two_point_synapses
      % error
      % Need to transform all annotations during window/region generation, pass in here...
      given_seg_raw_image = given_window_image(given_segmentation_image > 0);
      switch align_method
          case 'original'
            result_landmarks = [segmentation_2d_to_3d_landmark_function(given_segmentation_image, given_annotations(1, :)); ...
                               segmentation_2d_to_3d_landmark_function(given_segmentation_image, given_annotations(2, :))];
          case 'segment_centroid'
            result_landmarks = [segmentation_2d_to_3d_landmark_function_mean(given_segmentation_image, given_annotations(1, :)); ...
                               segmentation_2d_to_3d_landmark_function_mean(given_segmentation_image, given_annotations(2, :))];
          case 'segment_maximum'
            result_landmarks = [segmentation_2d_to_3d_landmark_function_max(given_segmentation_image, given_annotations(2, :)); ...
                               segmentation_2d_to_3d_landmark_function_max(given_segmentation_image, given_annotations(3, :))];
          case 'raw_original'
            result_landmarks = [segmentation_2d_to_3d_landmark_function(given_seg_raw_image, given_annotations(2, :)); ...
                               segmentation_2d_to_3d_landmark_function(given_seg_raw_image, given_annotations(3, :))];            
          case 'raw_centroid'
            result_landmarks = [segmentation_2d_to_3d_landmark_function_mean(given_seg_raw_image, given_annotations(2, :)); ...
                               segmentation_2d_to_3d_landmark_function_mean(given_seg_raw_image, given_annotations(3, :))];
          case 'raw_maximum'
            result_landmarks = [segmentation_2d_to_3d_landmark_function_max(given_seg_raw_image, given_annotations(2, :)); ...
                               segmentation_2d_to_3d_landmark_function_max(given_seg_raw_image, given_annotations(3, :))];
      end
      % keyboard
      % error
    else
      given_synapse_point = window_center_2d;
      result_landmarks = segmentation_2d_to_3d_landmark_function(given_segmentation_image, given_synapse_point);
    end
    
    centroid = image_centroid(given_segmentation_image);
    
    result_landmarks = [result_landmarks; centroid];
  end

  template_landmarks = [template_synapse; template_centroid];
  % warning('Adding two-point-related code!!!')
  if master_script_options.use_two_point_synapses
    % error
    template_landmarks = [template_synapse_segment; template_centroid];
    % % error
    % % Make the sides of the annotation 3D:
    % result_landmarks = [segmentation_2d_to_3d_landmark_function(given_segmentation_image, given_additional_annotations(1, :)); segmentation_2d_to_3d_landmark_function(given_segmentation_image, given_additional_annotations(2, :))];
  end
  
  debug_generate_figures = false;
  % debug_generate_figures = true
  % beep, keyboard
  % error('Implementation yet unfinished below this line!')
  
  t_cell_info.alignment_info = struct();
  t_cell_info.alignment_info.segmentation_landmark_function = @segmentation_landmark_function;
  t_cell_info.alignment_info.template_volume = template_volume;
  t_cell_info.alignment_info.segmentation_volume_threshold = segmentation_volume_threshold;
  %xruan 06/30/2015
  % t_cell_info.alignment_info.align_method = align_method;
  t_cell_info.alignment_info.align_method = master_script_options.align_method;
  
  
  if master_script_options.skip_alignment
    
    warning('>>>> HACK, options.skip_alignment true!')
    
  else    
    current_synapse_annotations = options.synapse_annotation;
    current_synapse_center_rounded = round([mean(current_synapse_annotations([1, 3])), mean(current_synapse_annotations([2, 4]))]);
    current_relative_time = relative_time;
    if isempty(current_synapse_center_rounded)
        fprintf('%s annotation is empty\n', current_synapse_annotations);
        return;
    end
    
    cell_alignment_filename = [sprintf('run%d_cell%02d_frame%05d_synapse%05d,%05d', frame_index(2), frame_index(3), frame_index(1), current_synapse_center_rounded)];;
    cell_segmentation_filename = [segmentations_filtered_location, cell_alignment_filename];
    cell_region_filename = [regions_location, cell_alignment_filename];
    [can_start, final_name, final_exists] = chunk_start_clean(alignments_location, cell_alignment_filename);
           
    if final_exists
        fprintf('Aligned cell %s already exists\n', cell_alignment_filename);
        return;
    end
    while ~final_exists
        if ~can_start && ~final_exists
            out_status = chunk_lock_clean(alignments_location);
            if out_status
                disp('Delete lock files whose jobs are not running!\n')
            end
        end
        try     
            a = load([cell_region_filename, '.mat']);
            current_synapse_annotations_transform = a.current_synapse_annotations_transform;            
            current_window_image = a.current_window_image;

            b = load([cell_segmentation_filename, '.mat']);
            segmentation_structure = b.segmentation_structure;
            window_number_slices = b.window_number_slices;
        catch
            save(final_name, '')
            chunk_finish(alignments_location, cell_alignment_filename);
            break;
        end

        current_synapse_annotations_transformed = reshape(current_synapse_annotations, 2, []).';
        current_synapse_annotations_transform_function = @(given_annotations)(given_annotations - repmat(current_synapse_annotations_transform.current_synapse_center_rounded, [size(given_annotations, 1), 1]) + current_synapse_annotations_transform.current_room) .* (current_synapse_annotations_transform.room / current_synapse_annotations_transform.current_room) + 1;
        % current_synapse_annotations_transform_function = @(given_annotations)(given_annotations - repmat(current_synapse_annotations_transform.current_synapse_center_rounded, [1, size(given_annotations, 2) / 2]) + current_synapse_annotations_transform.current_room) .* (current_synapse_annotations_transform.room / current_synapse_annotations_transform.current_room) + 1;
        current_synapse_annotations_transformed = current_synapse_annotations_transform_function(current_synapse_annotations_transformed);

        current_segmentation_image = segmentation_rasterization_function(segmentation_structure, [window_size_2d, window_number_slices]);
        % current_segmentation_image = segmentation_rasterization_function(segmentation_structure);
        % current_landmarks = segmentation_landmark_function(current_segmentation_image, current_synapse_annotations);
        align_method = master_script_options.align_method;
        current_landmarks = segmentation_landmark_function(current_window_image, current_segmentation_image, current_synapse_annotations_transformed, align_method);
        
        if isempty(current_landmarks)
          temp = struct; save(final_name, '-struct', 'temp'); chunk_finish(alignments_location, cell_alignment_filename);
          return;
        end
        % xruan 07/15/2015 check the relative end points of synapse to
        % centroid. The original coordinate system is left-handed, so
        % if the z coordinate of the the cross product is negative, it
        % means right end point is in the relative left, and need
        % change as left point. 

        if master_script_options.use_two_point_synapses
            CL = current_landmarks(1, :) - current_landmarks(3, :);
            CR = current_landmarks(2, :) - current_landmarks(3, :);
            cross_product = cross(CL, CR);
            if cross_product(3) < 0
                current_landmarks = current_landmarks([2, 1, 3], :);
            end
        end

        % Cell volume:
        current_volume = sum(current_segmentation_image(:));
        is_segmentation_low_volume = current_volume <= segmentation_volume_threshold;

        if any(~isfinite(current_landmarks(:)))
          % If the segmentation image is blank, current_landmarks will have nans:
          if debug_print_alignment_skip_reasons, fprintf('Skipping alignment %s current_synapse_center_rounded %d,%d frame %d, because no finite landmarks!\n', cell_alignment_filename, current_synapse_center_rounded, frame_index), end
           temp = struct; save(final_name, '-struct', 'temp'); chunk_finish(alignments_location, cell_alignment_filename);
          return
        end

        % At some point turn this off and filter in model building so can change this threshold to taste without recomputing:
        if is_segmentation_low_volume
          % If the segmentation image is blank, current_landmarks might not have nans?:
          if debug_print_alignment_skip_reasons, fprintf('Skipping alignment %s current_synapse_center_rounded %d,%d frame %d, because segmentation has low volume!\n', cell_alignment_filename, current_synapse_center_rounded, frame_index), end
          temp = struct; save(final_name, '-struct', 'temp'); chunk_finish(alignments_location, cell_alignment_filename);
          return
        end
        
        fprintf('Align segmented image for frame "%s"\n', cell_alignment_filename);

        % Compute absolute alignments between a frame and the template:
        try
          [~, ~, current_transform] = procrustes(template_landmarks, current_landmarks, 'Reflection', false, 'Scaling', false);
        catch raised_error
          getReport(raised_error, 'extended')
          % keyboard
          % rethrow(raised_error)
          temp = struct; save(final_name, '-struct', 'temp'); chunk_finish(alignments_location, cell_alignment_filename);
          return
        end
        current_transform_homogeneous = [diag(current_transform.b) * current_transform.T', current_transform.c(1, :)'; 0, 0, 0, 1];

        % Scale according to cell volume:
        current_volume = sum(current_segmentation_image(:));
        % Scale by volume about the synapse point:
        % current_transform_homogeneous = eye(4)
        % current_transform_homogeneous = translation_matrix3h(-current_landmarks(1, :)) * current_transform_homogeneous;
        % current_transform_homogeneous = translation_matrix3h(-current_landmarks(3, :)) * current_transform_homogeneous;
        current_transform_homogeneous = translation_matrix3h(-current_transform.c(1, :)') * current_transform_homogeneous;
        current_transform_homogeneous = scale_matrix3h(ones(1, 3) * nthroot(template_volume / current_volume, 3)) * current_transform_homogeneous;
        % current_transform_homogeneous = translation_matrix3h(current_landmarks(1, :)) * current_transform_homogeneous;
        % current_transform_homogeneous = translation_matrix3h(current_landmarks(3, :)) * current_transform_homogeneous;
        current_transform_homogeneous = translation_matrix3h(current_transform.c(1, :)') * current_transform_homogeneous;
        % current_transform_homogeneous = translation_matrix3h(current_landmarks(1, :)) * scale_matrix3h(ones(1, 3) * (template_volume / current_volume)) * translation_matrix3h(-current_landmarks(1, :)) * current_transform_homogeneous;

        current_transform_homogeneous_function = @(lambda, x_shift)[diag(current_transform.b) * current_transform.T' * lambda, ...
            current_transform.c(1, :)' + (1 - lambda) * diag(current_transform.b) * current_transform.T' * x_shift; ...
            0, 0, 0, 1];
        current_transform_homogeneous_old_function =  @(lambda, x_shift)[diag(current_transform.b) * current_transform.T' * lambda, ...
            lambda * current_transform.c(1, :)' + (1 - lambda) * x_shift; ...
            0, 0, 0, 1];

        current_volume = sum(current_segmentation_image(:));
        lambda = nthroot(template_volume / current_volume, 3);
        transformation_method_set = {'original', 'rescale_origin', 'rescale_centroid', 'rescale_left_end', 'rescale_right_end',...
            'old_left_end', 'old_right_end', 'old_centroid'};
        % transformation_method = 'original';
        % 02/22/2016 xruan
        % just use the method default to keep it consistent with the
        % running for the paper. 
        if master_script_options.use_two_point_synapses
            transformation_method = transformation_method_set{3};
        else
            % for one point synapse, still use centroid shift, but it
            % is current_landmarks(2, :);
            transformation_method = transformation_method_set{5};
        end

        switch transformation_method
            case 'original'
                current_transform_homogeneous = current_transform_homogeneous_function(1, [0, 0, 0]');
            case 'rescale_origin'
                current_transform_homogeneous = current_transform_homogeneous_function(lambda, [0, 0, 0]');
            case 'rescale_centroid'
                current_transform_homogeneous = current_transform_homogeneous_function(lambda, current_landmarks(3, :)');
            case 'rescale_left_end'
                current_transform_homogeneous = current_transform_homogeneous_function(lambda, current_landmarks(1, :)');
            case 'rescale_right_end'
                current_transform_homogeneous = current_transform_homogeneous_function(lambda, current_landmarks(2, :)');
            case 'old_left_end'
                current_transform_homogeneous = current_transform_homogeneous_old_function(lambda, current_landmarks(1, :)');
            case 'old_right_end'
                current_transform_homogeneous = current_transform_homogeneous_old_function(lambda, current_landmarks(2, :)');
            case 'old_centroid'
                current_transform_homogeneous = current_transform_homogeneous_old_function(lambda, current_landmarks(3, :)');
        end

        % Transform current_window_image:
        % xruan 06/30/2015
        % first pad the image to allow the transformed image to
        % maintain the signal
        % xruan 07/21/2015 
        % check if it is need to pad by compare the actual
        % window_number_slices with template slice numbers. if the
        % actual window number slices is larger, then pad 0. 
        padval = median(current_window_image(:));
        current_window_image_padded = padarray(current_window_image, [0, 0, max([0, cropped_size(3) - window_number_slices])], padval, 'post');
        current_window_image_transformed = affine(current_window_image_padded, current_transform_homogeneous([2, 1, 3, 4], [2, 1, 3, 4]), [], false, [], [], 'crop');
        % Transform segmentation_structure:
        segmentation_structure_transformed = segmentation_structure;
        segmentation_structure_transformed.vertices = current_transform_homogeneous * [segmentation_structure_transformed.vertices, ones(size(segmentation_structure_transformed.vertices, 1), 1)]';
        segmentation_structure_transformed.vertices = segmentation_structure_transformed.vertices(1:3, :)';

        % Crop the final aligned segmentation image:
        % cropped_segmentation_image_without_transformation = segmentation_rasterization_function(segmentation_structure, [window_size_2d, window_number_slices]);

        % xruan 06/30/2015
        %cropped_segmentation_image = segmentation_rasterization_function(segmentation_structure_transformed, [window_size_2d, window_number_slices]);
        cropped_segmentation_image = segmentation_rasterization_function(segmentation_structure_transformed, [window_size_2d, size(template_image, 3)]);
        cropped_segmentation_image_without_transformation = segmentation_rasterization_function(segmentation_structure, [window_size_2d,  size(template_image, 3)]);            % cropped_segmentation_image = segmentation_rasterization_function(segmentation_structure_transformed);

        cropped_segmentation_image_volume = sum(cropped_segmentation_image(:));

        cropped_segmentation_image = padarray(cropped_segmentation_image, max(cropped_size - size(cropped_segmentation_image), 0), 'post');
        cropped_segmentation_image = cropped_segmentation_image(1:cropped_size(1), 1:cropped_size(2), 1:cropped_size(3));

        % Set the edges to zero so the registration does not take forever:
        edge_width = 2;
        cropped_segmentation_image(1:edge_width, :, :) = 0;
        cropped_segmentation_image(end - edge_width + 1:end, :, :) = 0;
        cropped_segmentation_image(:, 1:edge_width, :) = 0;
        cropped_segmentation_image(:, end - edge_width + 1:end, :) = 0;
        cropped_segmentation_image(:, :, 1:edge_width) = 0;
        cropped_segmentation_image(:, :, end - edge_width + 1:end) = 0;

        % cropped_raw_image = current_window_image;
        cropped_raw_image = current_window_image_transformed;
        cropped_raw_image = padarray(cropped_raw_image, max(cropped_size - size(cropped_raw_image), 0), 'post');
        cropped_raw_image = cropped_raw_image(1:cropped_size(1), 1:cropped_size(2), 1:cropped_size(3));

        % if debug_keyboard_alignment_skips && current_volume ~= cropped_segmentation_image_volume
          % dbstack, beep, keyboard
        % end

        is_cropped_segmentation_low_volume = cropped_segmentation_image_volume <= segmentation_volume_threshold;
        if is_cropped_segmentation_low_volume
          if debug_print_alignment_skip_reasons, disp(sprintf('Skipping alignment %s current_synapse_center_rounded %d,%d frame %d, because segmentation has low volume!\n', run_key, current_synapse_center_rounded, frame_index)), end
          temp = struct; save(final_name, '-struct', 'temp'); chunk_finish(alignments_location, cell_alignment_filename);
          return
        end

        % At some point save current_cell_volumes:
        save(final_name, 'current_landmarks', 'current_transform', 'current_transform_homogeneous', 'cropped_segmentation_image', 'cropped_raw_image')
        chunk_finish(alignments_location, cell_alignment_filename);
        final_exists = true;
     end  
  end
end

