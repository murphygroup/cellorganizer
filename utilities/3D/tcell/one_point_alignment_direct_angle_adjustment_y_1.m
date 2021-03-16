function [theta_y_adjust] = one_point_alignment_direct_angle_adjustment_y_1(current_landmarks, landmark_vector, Thetas_2pt, euler_theta_y, infer_model_param)

%
% modularize the process for adjustment of rotation angle by y-axis. 
%
% idea: learn a Gaussian mixture model of y angles, and then fit to the model. 
% 
%
%
% Author: Xiongtao Ruan
% Date: Sept. 13, 2016

k = 11;

current_landmarks_vector = current_landmarks(1, :) - current_landmarks(2, :);
 
current_landmarks_vector = current_landmarks_vector ./ norm(current_landmarks_vector);

landmark_vector = bsxfun(@rdivide, landmark_vector, sqrt(sum(landmark_vector .^ 2, 2))); 

vector_distance = sqrt(sum(bsxfun(@minus, landmark_vector, current_landmarks_vector) .^ 2, 2));

[~, sorted_inds] = sort(vector_distance);

% voting

k_sorted_inds = sorted_inds(1  : k);

k_sorted_dists = vector_distance(k_sorted_inds);
k_Thetas_2pt_y = Thetas_2pt(k_sorted_inds, 2);

k_labels = zeros(k, 1);

k_labels = (abs(k_Thetas_2pt_y) > abs(k_Thetas_2pt_y + pi)) + (abs(k_Thetas_2pt_y) >= abs(k_Thetas_2pt_y - pi));

cur_label = 0;

if sum(k_labels) > k / 2
    cur_label = 1;
end

k_Thetas_2pt_yl = k_Thetas_2pt_y(k_labels == cur_label);
k_sorted_dists_l = k_sorted_dists(k_labels == cur_label);

% angles around -pi and pi would mess up each other

if cur_label == 1
    k_Thetas_2pt_yl(k_Thetas_2pt_yl < 0) = k_Thetas_2pt_yl(k_Thetas_2pt_yl < 0) + 2 * pi;
end    

theta_y_average = (1 - k_sorted_dists_l)' * k_Thetas_2pt_yl / sum(1 - k_sorted_dists_l);

theta_y_adjust = theta_y_average - euler_theta_y;

if theta_y_adjust < -pi
    theta_y_adjust = theta_y_adjust + 2 * pi;
end

if theta_y_adjust > 1.5 * pi
    theta_y_adjust = theta_y_adjust - 2 * pi;
end

end