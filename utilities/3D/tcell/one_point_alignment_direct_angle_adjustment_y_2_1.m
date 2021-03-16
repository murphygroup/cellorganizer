function [theta_y_adjust] = one_point_alignment_direct_angle_adjustment_y_2_1(I_seg, I_raw, synapse_center,  cur_landmark_vector, Thetas_2pt, euler_theta_y, landmark_vector,  Btree)
%
% idea: modularize the process for adjustment of rotation angle by y-axis. 
%
%
% Author: Xiongtao Ruan
% Date: Sept. 14, 2016


% decide whether rotation pi/2 by y-axis
% [I_raw_transformed_z, I_seg_transformed_z] = image_transformation_from_euler_angles(I_raw, I_seg, synapse_center, 0, -pi / 2, 0);

if abs(sin(euler_theta_y)) >= sqrt(2) / 2
    [I_raw_transformed_z, I_seg_transformed_z] = image_transformation_from_euler_angles(I_raw, I_seg, synapse_center, 0, -pi / 2, 0);
    I_seg_transformed = I_seg_transformed_z;
    I_raw_transformed = I_raw_transformed_z;
else
    I_seg_transformed = I_seg;
    I_raw_transformed = I_raw;    
end


I_raw_cropped = I_raw_transformed .* (I_seg_transformed >= 0.5);
sz = size(I_raw_cropped, 3);
fs_1 = sum(reshape(I_raw_cropped, [], sz));
fs_3 = [0, 0];
fs_3(1) = sum(reshape(I_raw_cropped(:, :, 1 : floor(sz / 2)), 1, []));
fs_3(2) = sum(reshape(I_raw_cropped(:, :, ceil(sz / 2) + 1 : end), 1, []));


cur_landmark_vector = cur_landmark_vector / norm(cur_landmark_vector);
landmark_vector = bsxfun(@rdivide, landmark_vector, sqrt(sum(landmark_vector .^ 2, 2)));
Thetas_2pt_y = Thetas_2pt(:, 2);
Theta_2y_labels = (abs(Thetas_2pt_y) > abs(Thetas_2pt_y + pi)) + (abs(Thetas_2pt_y) >= abs(Thetas_2pt_y - pi));

Xte = [fs_1 / max(fs_1, [], 2),  fs_3 / max(fs_3, [], 2), cur_landmark_vector ./ norm(cur_landmark_vector)];

Yfit = predict(Btree, Xte);
cur_label = cellfun(@str2num, Yfit);

theta_y_adjust = 0;
if cur_label == 1
    theta_y_adjust = pi;
end

theta_y_adjust = theta_y_adjust - euler_theta_y;

if theta_y_adjust > 1.25 * pi
    theta_y_adjust = theta_y_adjust - 2 * pi;
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
