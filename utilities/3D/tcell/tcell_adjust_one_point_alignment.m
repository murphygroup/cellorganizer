function [t_cell_info] = tcell_adjust_one_point_alignment(t_cell_info, options)
% the function is used for the adjustment of one point alignment based on
% model of two point alignment.
% 
% Copied form master_script_adjust_one_point_alignment.m
%
% Author: Xiongtao Ruan
% Date: 08/21/2018



% Variables referenced more than once can be copied to local variables:
master_script_options = t_cell_info.options;
regions_location = t_cell_info.path_info.regions_location;
segmentations_filtered_location = t_cell_info.path_info.segmentations_filtered_location;
alignments_location = t_cell_info.path_info.alignments_location;
alignments_adjust_location = t_cell_info.path_info.alignments_adjust_location; 

template_centroid = t_cell_info.template_info.template_centroid;
template_synapse = t_cell_info.template_info.template_synapse;
template_image = t_cell_info.template_info.template_image;
number_slices = t_cell_info.template_info.template_number_slices;

window_size_2d = t_cell_info.preprocessing_info.window_size_2d;
window_center_2d = t_cell_info.preprocessing_info.window_center_2d;
cropped_size = t_cell_info.preprocessing_info.cropped_size;

segmentation_rasterization_function = t_cell_info.segmentation_info.segmentation_rasterization_function;

template_volume = t_cell_info.alignment_info.template_volume;
segmentation_volume_threshold = t_cell_info.alignment_info.segmentation_volume_threshold;

one_point_alignment_adjust_model_filename = t_cell_info.options.one_point_alignment_adjust_model_filename;


image_name = options.image_name;
synapse_annotation = options.synapse_annotation;
relative_time = options.relative_time;
frame_channel = options.frame_channel;
frame_index = options.frame_index;
run_index = frame_index(:, 2);  

% debug_print_alignments_adjust_skip_reasons = false;
debug_print_alignments_adjust_skip_reasons = true;

% check_morph_file_loadility = master_script_options.check_morph_file_loadility;
% Set options for T cell image registration:

function [raw_image_transformed, seg_image_transformed, T] = image_transformation_from_euler_angles(raw_image, seg_image, theta_x, theta_y, theta_z, angle_order)
    if nargin < 6
        angle_order = 'ZXY';
    end
    [M] = euler_angles_to_tranformation_matrix(theta_x, theta_y, theta_z, 1, angle_order);
    T = eye(4);
    T(1 : 3, 1 : 3) = M;
    % rotate using template centroid as origin
    T(1 : 3, 4) = (eye(3) - M) * template_synapse';            
    seg_image_transformed = affine(seg_image, T([2, 1, 3, 4], [2, 1, 3, 4]), [], false, [], [], 'crop');
    raw_image_transformed = affine(raw_image, T([2, 1, 3, 4], [2, 1, 3, 4]), [], false, [], [], 'crop');
end

% load regression and classification model
persistent beta beta_y mvr_param Btree euler_angle_order infer_model_param ...
    Thetas_1pt Thetas_2pt landmark_vector_1pt landmark_vector_2pt raw_landmark_vector_1pt raw_landmark_vector_2pt

if isempty(beta)
    b = load(one_point_alignment_adjust_model_filename);
    beta = b.model1.beta;
    beta_y = b.model3.beta;
    mvr_param = b.model1.mvr_param;

    Btree = b.model2.Btree;
    % Btree_2 = b.model3.Btree;

    % rf_param = b.model2.param;
    euler_angle_order = b.param.angle_order;
    infer_model_param = b.param;
    Thetas_1pt = b.model_info.Thetas_1pt;
    Thetas_2pt = b.model_info.Thetas_2pt;
    landmark_vector_1pt = b.model_info.landmark_vector_1pt;
    landmark_vector_2pt = b.model_info.landmark_vector_2pt;
    raw_landmark_vector_1pt = b.model_info.raw_landmark_vector_1pt;
    raw_landmark_vector_2pt = b.model_info.raw_landmark_vector_2pt;
end

% Keep running until everything is finished:
current_synapse_annotations = options.synapse_annotation;
current_synapse_center_rounded = round([mean(current_synapse_annotations([1, 3])), mean(current_synapse_annotations([2, 4]))]);
current_synapse_center_rounded = current_synapse_center_rounded(1, 1:2);
current_relative_time = relative_time;
if isempty(current_synapse_center_rounded)
    fprintf('%s annotation is empty\n', current_synapse_annotations);
    return;
end

cell_alignment_adjust_filename = [sprintf('run%d_cell%02d_frame%05d_synapse%05d,%05d', frame_index(2), frame_index(3), frame_index(1), current_synapse_center_rounded)];
cell_segmentation_filename = [segmentations_filtered_location, cell_alignment_adjust_filename];
cell_region_filename = [regions_location, cell_alignment_adjust_filename];
cell_alignment_filename = [alignments_location, cell_alignment_adjust_filename];
[can_start, final_name, final_exists] = chunk_start_clean(alignments_adjust_location, cell_alignment_adjust_filename);

if ~exist([cell_alignment_filename, '.mat'], 'file')
  % Given modifications to master_script_align_segmentations.m, all alignment should already be complete and thus this file will never be available (should look into why it needs to be checked, though):
end

% frame_count = 0;
% filenames = {};
% Theta_y_infer = [];
% Xt = [];
% Theta_y_infer_final = [];

while ~final_exists
    if ~can_start && ~final_exists
        out_status = chunk_lock_clean(alignments_location);
        if out_status
            disp('Delete lock files whose jobs are not running!\n')
        end
    end


    if final_exists && false
        a = load(final_name)
        % Debug plot:
        % imshow(reshape_2d([contrast_stretch(a.cropped_raw_image); a.cropped_segmentation_image; template_image; a.cropped_segmentation_image_deformed], -1)), pause
        % keyboard
        % imshow(reshape_2d([template_image; contrast_stretch(a.cropped_raw_image); contrast_stretch(a.cropped_raw_image_deformed); a.cropped_segmentation_image; a.cropped_segmentation_image_deformed], -1)), pause
        imshow(reshape_2d([template_image; a.cropped_segmentation_image; a.cropped_segmentation_image_deformed; contrast_stretch(a.cropped_raw_image); contrast_stretch(a.cropped_raw_image_deformed)], -1)), pause
        % imshow(reshape_2d([template_image; a.cropped_segmentation_image; a.cropped_segmentation_image_deformed; contrast_stretch(a.cropped_raw_image); contrast_stretch(a.cropped_raw_image_deformed)], -1)), pause
        % imshow(reshape_contrast([a.cropped_raw_image; a.cropped_raw_image .* a.cropped_segmentation_image], -1)), pause
        % imshow(reshape_contrast([a.cropped_raw_image_deformed; a.cropped_raw_image_deformed .* a.cropped_segmentation_image_deformed], -1)), pause
        % imshow(reshape_contrast(a.cropped_raw_image_deformed, -1)), pause
        % imshow(reshape_2d([; contrast_stretch([a.cropped_raw_image_deformed; a.cropped_raw_image_deformed .* a.cropped_segmentation_image_deformed])], -1)), pause
    end

    if ~can_start
        continue
    end

    a = load([cell_segmentation_filename, '.mat']);
    % Feb. 14, 2016 xruan disable it
    % segmentation_structure = a.segmentation_structure;
    % try
    a = load([cell_alignment_filename, '.mat']);
    % catch raised_error
    % getReport(raised_error, 'extended')
    % end
    if isempty(fieldnames(a))
        temp = struct; save(final_name, '-struct', 'temp'); chunk_finish(alignments_adjust_location, cell_alignment_adjust_filename);
        fprintf('Skipping alignment adjustment for frame "%s" because alignment empty\n', cell_alignment_adjust_filename);
        return
    end
    current_landmarks = a.current_landmarks;
    current_transform = a.current_transform;
    % xruan 08/14/2015
    % for wave2 active data, make sure it is double for deformation. 
    cropped_segmentation_image = double(a.cropped_segmentation_image);  
    % cropped_segmentation_image = a.cropped_segmentation_image;
    cropped_raw_image = a.cropped_raw_image;
    current_transform_homogeneous = a.current_transform_homogeneous;

    current_volume = sum(cropped_segmentation_image(:));
    is_segmentation_low_volume = current_volume <= segmentation_volume_threshold;

    if is_segmentation_low_volume
    % if debug_print_alignments_adjust_skip_reasons, fprintf_diary_cycle_function('Skipping morph %s cell %d because cropped_segmentation_image blank\n', run_key, current_cell_index), end
        fprintf('Skipping alignment adjustment for frame "%s" because of bad segmentation\n', cell_alignment_adjust_filename);
        % chunk_finish(deformations_location, cell_alignment_adjust_filename);
        temp = struct; save(final_name, '-struct', 'temp'); chunk_finish(alignments_adjust_location, cell_alignment_adjust_filename);
        return
    end

    fprintf('Adjust aligned segmentations for frame "%s"\n', cell_alignment_adjust_filename)

    if false && ~isempty(regexp(cell_alignment_adjust_filename, 'Actin Regulator - Full Stimulus_Actin_Actin 02 24 14 run 1.csv_frame00008_synapse00279,00045', 'once'))
        flag = 1;
    end                

    [euler_theta_x, euler_theta_y, euler_theta_z] = tranformation_matrix_to_euler_angles(current_transform_homogeneous, euler_angle_order);
    feature_param = struct();   
    [cur_feature] = one_point_alignment_feature_extraction(a, feature_param);           
    cur_feature = [cur_feature', 1];
    % infere theta_x, and theta_z
    Y_rec = cur_feature * beta;

    if strcmp(mvr_param.method, 'SinCos_ZX')
        cur_param = [Y_rec(1), 0, Y_rec(2), Y_rec(3), 1, Y_rec(4)];
    end

    theta_1 = atan2(cur_param(1), cur_param(4));
    theta_2 = atan2(cur_param(2), cur_param(5));
    theta_3 = atan2(cur_param(3), cur_param(6));

    [cropped_raw_image_zx, cropped_segmentation_image_zx, current_adjust_transform_zx] = image_transformation_from_euler_angles(cropped_raw_image, cropped_segmentation_image, theta_1, 0, theta_3);                
    [~, euler_theta_y_zx] = tranformation_matrix_to_euler_angles(current_adjust_transform_zx * current_transform_homogeneous, euler_angle_order);

    % options.infer_theta_y_method = 'method_1';
    switch master_script_options.infer_theta_y_method
        case 'method_0'
            [theta_y_adjust] = one_point_alignment_direct_angle_adjustment_y_0(cropped_segmentation_image_zx, cropped_raw_image_zx, template_synapse, euler_theta_y_zx, infer_model_param)
        case 'method_1'
            [theta_y_adjust] = one_point_alignment_direct_angle_adjustment_y_1(current_landmarks, landmark_vector_1pt, Thetas_2pt, euler_theta_y_zx, infer_model_param);
        case 'method_2'
            [theta_y_adjust] = one_point_alignment_direct_angle_adjustment_y_1(current_landmarks, landmark_vector_2pt, Thetas_2pt, euler_theta_y_zx, infer_model_param);
        case 'method_3'
            current_landmarks_vector = current_landmarks(1, :) - current_landmarks(2, :);
            [theta_y_adjust, Xte] = one_point_alignment_direct_angle_adjustment_y_2(cropped_segmentation_image, cropped_raw_image, current_landmarks_vector, Thetas_2pt, euler_theta_y, Btree);
            [cropped_raw_image_y, cropped_segmentation_image_y, current_adjust_transform] = image_transformation_from_euler_angles(cropped_raw_image_zx, cropped_segmentation_image_zx, 0, theta_y_adjust, 0);
            current_adjust_transform = current_adjust_transform * current_adjust_transform_zx;
        case 'method_4'
            cur_param_y = cur_feature * beta_y;
            theta_2 = atan2(cur_param_y(1), cur_param_y(2));
            theta_y_adjust = theta_2 - euler_theta_y;
            [cropped_raw_image_y, cropped_segmentation_image_y, current_adjust_transform] = image_transformation_from_euler_angles(cropped_raw_image_zx, cropped_segmentation_image_zx, 0, theta_y_adjust, 0);
            current_adjust_transform = current_adjust_transform * current_adjust_transform_zx;
        case 'method_5'
            current_landmarks_vector = current_landmarks(1, :) - current_landmarks(2, :);
            [theta_y_adjust] = one_point_alignment_direct_angle_adjustment_y_2_2(cropped_segmentation_image, cropped_raw_image, template_synapse, current_landmarks_vector, Thetas_2pt, euler_theta_y, landmark_vector_1pt, Btree);
            [cropped_raw_image_y, cropped_segmentation_image_y, current_adjust_transform] = image_transformation_from_euler_angles(cropped_raw_image_zx, cropped_segmentation_image_zx, 0, theta_y_adjust, 0);
            current_adjust_transform = current_adjust_transform * current_adjust_transform_zx;                        
    end

    if false
        sz = size(cropped_segmentation_image_y, 3);
        fs_9_1 = sum(reshape(cropped_segmentation_image_y(:, :, 1 : floor(sz / 2)), 1, []));
        fs_9_2 = sum(reshape(cropped_segmentation_image_y(:, :, ceil(sz / 2) + 1 : end), 1, []));

        if fs_9_2 / fs_9_1 > 2 
            figure, imshow(reshape_2d(cropped_raw_image_y), []);
            figure, imshow(reshape_2d(cropped_segmentation_image_y), []);
            flag =  1
        end
    end

    if master_script_options.infer_theta_z_separately 

        [opt_angles] = one_point_alignment_direct_angle_adjustment_1(cropped_segmentation_image_y, cropped_raw_image_y, template_synapse, infer_model_param);

        [cropped_raw_image_z, cropped_segmentation_image_z, current_adjust_transform_1] = image_transformation_from_euler_angles(cropped_raw_image_y, cropped_segmentation_image_y, opt_angles(1), 0, opt_angles(2));
    else
        current_adjust_transform_1 = eye(4);
        cropped_raw_image_z = cropped_raw_image_y;
        cropped_segmentation_image_z = cropped_segmentation_image_y;                                        
    end

    if false
        [opt_angles_1] = one_point_alignment_direct_angle_adjustment_2(cropped_segmentation_image_z, cropped_raw_image_z, template_synapse, infer_model_param);

        [cropped_raw_image, cropped_segmentation_image, current_adjust_transform_2] = image_transformation_from_euler_angles(cropped_raw_image_z, cropped_segmentation_image_z, opt_angles_1(1), opt_angles(2), 0);
    else
        current_adjust_transform_2 = eye(4);
        cropped_raw_image = cropped_raw_image_z;
        cropped_segmentation_image = cropped_segmentation_image_z;

    end

    current_adjust_transform = current_adjust_transform_2 * current_adjust_transform_1 * current_adjust_transform;
    % decide whether flip the image after deciding rotation, that is
    % whether rotate another pi by y-axis

    if false
        frame_count = frame_count + 1;
        filenames = [filenames; {cell_alignment_adjust_filename}];
        Theta_y_infer = [Theta_y_infer; theta_y_adjust + euler_theta_y_zx]; 
        [~, theta_2_final] = tranformation_matrix_to_euler_angles(current_adjust_transform * current_transform_homogeneous, euler_angle_order);
        Theta_y_infer_final = [Theta_y_infer_final; theta_2_final];
        Xt = [Xt; Xte];
        if abs(abs(theta_y_adjust + euler_theta_y_zx) - abs(theta_2_final)) > 0.5
            flag = 1;
        end
    end   

    if false
    % Debug plot:
        figure,
        image_to_show = [];
        image_to_show = [image_to_show; cropped_raw_image_y];
        image_to_show = [image_to_show; cropped_raw_image];              
        imshow(reshape_2d(image_to_show), []);
    end

    save(final_name, 'current_landmarks', 'current_transform', 'current_transform_homogeneous', 'current_adjust_transform', 'cropped_segmentation_image', 'cropped_raw_image')
    chunk_finish(alignments_adjust_location, cell_alignment_adjust_filename);
    final_exists = true;    
end

end


