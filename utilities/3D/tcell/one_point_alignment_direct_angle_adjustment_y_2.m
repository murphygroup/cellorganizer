function [theta_y_adjust, Xte] = one_point_alignment_direct_angle_adjustment_y_2(I_seg, I_raw, landmark_vector, Thetas_2pt, euler_theta_y, Btree)
%
% idea: modularize the process for adjustment of rotation angle by y-axis. 
%
%
% Author: Xiongtao Ruan
% Date: Sept. 14, 2016

I_raw_cropped = I_raw .* (I_seg >= 0.5);
sz = size(I_raw_cropped, 3);
fs_1 = sum(reshape(I_raw_cropped, [], sz));
fs_3 = [0, 0];
fs_3(1) = sum(reshape(I_raw_cropped(:, :, 1 : floor(sz / 2)), 1, []));
fs_3(2) = sum(reshape(I_raw_cropped(:, :, ceil(sz / 2) + 1 : end), 1, []));

Xte = [fs_1 / max(fs_1, [], 2),  fs_3 / max(fs_3, [], 2), landmark_vector ./ norm(landmark_vector)];

Yfit = predict(Btree, Xte);
flip_label = cellfun(@str2num, Yfit);

theta_y_adjust = 0;
if flip_label == 1
    theta_y_adjust = pi;
end

theta_y_adjust = theta_y_adjust - euler_theta_y;

if theta_y_adjust > 1.5 * pi
    theta_y_adjust = theta_y_adjust - 2 * pi;
end

end

