function [landmarks_shape, scales, centers, angles, cell_nuc_scales, cell_nuc_centers, cell_nuc_angles] = pca_shape_processing(all_landmarks_cell, options)
% the overall function for processing cell and nuclear landmarks. 
% 
% Author: Xiongtao Ruan
% Date: 01/13/2019
%
%


if numel(all_landmarks_cell) == 1
    pca_model_type = 'separate';
elseif numel(all_landmarks_cell) == 2
    pca_model_type = 'joint';
end
    

options = ml_initparam(options, struct('shape_type', 'shape_size', ...  % other choices: shape, preshape, original
                                       'alignment_method', 'cell_align', ...  % other choices: separate_align, joint align
                                       'outline_scaling_method', 'joint_scaling')); % other choices: separate_scaling

if strcmp(pca_model_type, 'separate')
    cell_landmarks = all_landmarks_cell{1};
    switch options.shape_type
        case 'original'
            landmarks_shape = cell_landmarks;
            scales = [];
            centers = [];
            angles = [];
        case 'preshape'
            [landmarks_shape, scales, centers] = convert_to_preshape(landmarks);
            angles = [];
        case {'shape', 'shape_size'}
            [landmarks_shape, scales, centers, angles] = convert_to_shape(cell_landmarks);
    end
elseif strcmp(pca_model_type, 'joint')
    cell_landmarks = all_landmarks_cell{1};
    nuc_landmarks = all_landmarks_cell{2};

    switch options.alignment_method
        case 'separate_align'
            [cell_landmarks_shape, cell_scales, cell_centers, cell_angles] = convert_to_shape(cell_landmarks);
            [nuc_landmarks_shape, nuc_scales, nuc_centers, nuc_angles] = convert_to_shape(nuc_landmarks);
            landmarks_shape = cat(1, cell_landmarks_shape, nuc_landmarks_shape);
            scales = cell_scales;
            centers = cell_centers;
            angles = cell_angles;
        case 'cell_align'
            [cell_landmarks_shape, cell_scales, cell_centers, cell_angles] = convert_to_shape(cell_landmarks);

            cur_options = struct();
            cur_options.use_given_center = true;
            cur_options.centers = cell_centers;
            cur_options.use_given_rotation_matrix = true;
            cur_options.angles = cell_angles;
            [nuc_landmarks_shape, nuc_scales, nuc_centers, nuc_angles] = convert_to_shape(nuc_landmarks, cur_options);

            landmarks_shape = cat(1, cell_landmarks_shape .* cell_scales, nuc_landmarks_shape .* nuc_scales);
            if strcmp(options.outline_scaling_method, 'joint_scaling')
                scales = sqrt(sum(sum(landmarks_shape .^ 2, 1), 2));
                landmarks_shape = landmarks_shape ./ scales;
                cell_landmarks_shape = cell_landmarks_shape .* cell_scales ./ scales;
                nuc_landmarks_shape = nuc_landmarks_shape .* nuc_scales ./ scales;

                cell_scales = scales;
                nuc_scales = scales;
                centers = cell_centers;
                angles = cell_angles;
            else
                scales = cat(4, cell_scales, nuc_scales);
                centers = cell_centers;
                angles = cell_angles;
            end
        case 'joint_align'
            N_cell_points = size(all_landmarks_cell{1}, 1);
            all_landmarks = cat(1, cell_landmarks, nuc_landmarks);
            cur_options = struct();
            cur_options.N_cell_points = N_cell_points;
            cur_options.component = 'both';
            [landmarks_shape, scales, centers, angles] = convert_to_shape(all_landmarks, cur_options);
            cell_landmarks_shape = landmarks_shape(1 : N_cell_points, :, :);
            nuc_landmarks_shape = landmarks_shape(N_cell_points + 1 : end, :, :);
            cell_scales = scales;
            cell_centers = centers;
            cell_angles = angles;
            nuc_scales = scales;
            nuc_centers = centers;
            nuc_angles = angles;
    end
    cell_nuc_scales = {cell_scales, nuc_scales};
    cell_nuc_centers = {cell_centers, nuc_centers};
    cell_nuc_angles = {cell_angles, nuc_angles};
end


end