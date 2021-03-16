function [theta_y_adjust] = one_point_alignment_direct_angle_adjustment_y_2_2(I_seg, I_raw, synapse_center,  cur_landmark_vector, Thetas_2pt, euler_theta_y, landmark_vector,  Btree)
%
% idea: modularize the process for adjustment of rotation angle by y-axis. 
%
%
% Author: Xiongtao Ruan
% Date: Sept. 14, 2016

cur_landmark_vector = cur_landmark_vector / norm(cur_landmark_vector);
landmark_vector = bsxfun(@rdivide, landmark_vector, sqrt(sum(landmark_vector .^ 2, 2)));
Thetas_2pt_y = Thetas_2pt(:, 2);
Theta_2y_labels = (abs(Thetas_2pt_y) > abs(Thetas_2pt_y + pi)) + (abs(Thetas_2pt_y) >= abs(Thetas_2pt_y - pi));

cur_label = 0;

landmark_vector_2pt_label = landmark_vector(Theta_2y_labels == cur_label, :);
Thetas_2pt_y_label = Thetas_2pt_y(Theta_2y_labels == cur_label);
k = 5;
vector_distance = sqrt(sum(bsxfun(@minus, landmark_vector_2pt_label, cur_landmark_vector) .^ 2, 2));
[~, sorted_inds] = sort(vector_distance);
k_sorted_inds = sorted_inds(1  : k);
k_sorted_dists = vector_distance(k_sorted_inds);
k_Thetas_2pt_y = Thetas_2pt_y_label(k_sorted_inds);
k_labels = (abs(k_Thetas_2pt_y) > abs(k_Thetas_2pt_y + pi)) + (abs(k_Thetas_2pt_y) >= abs(k_Thetas_2pt_y - pi));
k_Thetas_2pt_yl = k_Thetas_2pt_y(k_labels == cur_label);
k_sorted_dists_l = k_sorted_dists(k_labels == cur_label);

if cur_label == 1
    k_Thetas_2pt_yl(k_Thetas_2pt_yl < 0) = k_Thetas_2pt_yl(k_Thetas_2pt_yl < 0) + 2 * pi;
end    

% change to use Gaussian weight
sigma = median(k_sorted_dists_l);
theta_y_average = exp(- 1 / 2 * k_sorted_dists_l .^2 / sigma^2)' * k_Thetas_2pt_yl / sum(exp(- 1/ 2 * k_sorted_dists_l .^ 2/ sigma^2));

theta_y_adjust = theta_y_average - euler_theta_y;

if theta_y_adjust < -pi
    theta_y_adjust = theta_y_adjust + 2 * pi;
end

if theta_y_adjust > 1.25 * pi
    theta_y_adjust = theta_y_adjust - 2 * pi;
end

end
