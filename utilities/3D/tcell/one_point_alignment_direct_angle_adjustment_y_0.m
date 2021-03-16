function [opt_angles] = one_point_alignment_direct_angle_adjustment_y_0(I_seg, I_raw, synapse_center, euler_theta_y, param)
%
% idea: modularize the process for adjustment of rotation angle by y-axis. 
%
%
% Author: Xiongtao Ruan
% Date: Sept. 13, 2016

[center_slice_roundness, rms] = rotation_by_y_axis_classification_feature_extration(I_raw, I_seg);

% decide whether rotation pi/2 by y-axis
[I_raw_transformed_z, I_seg_transformed_z] = image_transformation_from_euler_angles(I_raw, I_seg, synapse_center, 0, -pi / 2, 0);

[center_slice_roundness_z, rms_z] = rotation_by_y_axis_classification_feature_extration(I_raw_transformed_z, I_seg_transformed_z);

% use random forest model to predict whether rotate or not. 
X_rotation = [center_slice_roundness, center_slice_roundness_z, rms, rms_z, euler_theta_y];
Yfit = predict(Btree, X_rotation);
rotation_label = cellfun(@str2num, Yfit);

theta_y_adjust = 0;
if rotation_label == 1
    I_seg_transformed = I_seg_transformed_z;
    I_raw_transformed = I_raw_transformed_z;
    theta_y_adjust = -pi / 2;
else
    I_seg_transformed = I_seg;
    I_raw_transformed = I_raw;    
end

% compute the percentile of each slice to decide whether flip the
% image
slices_prctile_thresh = 90;        
for k = 1 : floor((number_slices + 1) ./ 2)
    cur_prctile = zeros(1, 2);
    cur_raw_slice = squeeze(I_raw_transformed(:, :, k));
    cur_seg_slice = squeeze(I_seg_transformed(:, :, k));
    cur_raw_slice_1 = squeeze(I_raw_transformed(:, :, number_slices - k + 1));
    cur_seg_slice_1 = squeeze(I_seg_transformed(:, :, number_slices - k + 1));                

    if ~isempty(cur_raw_slice(cur_seg_slice >= 0.5))
        cur_prctile(1) = prctile(cur_raw_slice(cur_seg_slice >= 0.5), slices_prctile_thresh);
    end              

    if ~isempty(cur_raw_slice_1(cur_seg_slice_1 >= 0.5))
        cur_prctile(2) = prctile(cur_raw_slice_1(cur_seg_slice_1 >= 0.5), slices_prctile_thresh);
    end

    percentiles(k, :) = cur_prctile;
end
percentiles(percentiles(:, 1) == 0, :) = [] ;
prticle_ratio = sum(percentiles(:, 1) >= percentiles(:, 2)) ./size(percentiles, 1);

prctile_ratio_thresh = 0.35;
if prticle_ratio <= prctile_ratio_thresh;
    theta_y_adjust = theta_y_adjust + pi;
end

end


function [center_slice_roundness, rms] = rotation_by_y_axis_classification_feature_extration(cropped_raw_image_transformed, cropped_segmentation_image_transformed)

    center_slice_number = floor((size(cropped_raw_image_transformed, 3) + 1) / 2);
    I_raw = squeeze(cropped_raw_image_transformed(:, :, center_slice_number) .* (cropped_segmentation_image_transformed(:, :, center_slice_number) >= 0.5)); 
    I_1 = I_raw .* (I_raw == 0) .* min(I_raw(I_raw > 0)) + I_raw;
    I_2 = (I_1 - min(I_1(:))) ./ (max(I_1(:)) - min(I_1(:)));
    thresh = graythresh(I_2);
    I_2_seg = I_2 >= thresh;
    states = regionprops(I_2_seg, 'MajorAxisLength', 'MinorAxisLength', 'Orientation');
    try
        center_slice_roundness = states.MajorAxisLength /  states.MinorAxisLength;
    catch
        center_slice_roundness = 1;
    end

    two_side_slice_number = 8;
    z_range = center_slice_number + (-two_side_slice_number : two_side_slice_number);
    rms = [var(reshape(cropped_raw_image_transformed(:, :, z_range) .*(cropped_segmentation_image_transformed(:, :, z_range) >= 0.5), 1, []))];

end

function [raw_image_transformed, seg_image_transformed, T] = image_transformation_from_euler_angles(raw_image, seg_image, rotation_center, theta_x, theta_y, theta_z, angle_order)
    if nargin < 7
        angle_order = 'ZXY';
    end
    [M] = euler_angles_to_tranformation_matrix(theta_x, theta_y, theta_z, 1, angle_order);
    T = eye(4);
    T(1 : 3, 1 : 3) = M;
    % rotate using template centroid as origin
    T(1 : 3, 4) = (eye(3) - M) * rotation_center';            
    seg_image_transformed = affine(seg_image, T([2, 1, 3, 4], [2, 1, 3, 4]), [], false, [], [], 'crop');
    raw_image_transformed = affine(raw_image, T([2, 1, 3, 4], [2, 1, 3, 4]), [], false, [], [], 'crop');
end
