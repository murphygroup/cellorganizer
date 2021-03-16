function [landmarks_out, scales, centers, angles] = convert_to_shape(landmarks, options)
% convert the landmarks to shape and the same time save the scale factor and centers. 
% 05/25/2018 also make the starting point the same location, i.e. with
% smallest x-axis

% landmarks_c = squeeze(landmarks(:, 1, :) + 1i * landmarks(:, 2, :));

if nargin < 2
    options = struct();
end

options = process_options_structure(struct('component', 'cell', ...
                                           'N_cell_points', 2000), options);

component = options.component;
if strcmp(component, 'both')
    N_cell_points = options.N_cell_points;
end

if isfield(options, 'use_given_center') && options.use_given_center
    centers = options.centers;
else
	centers = mean(landmarks, 1);
end
landmarks_c1 = landmarks - centers;

scales = sqrt(sum(sum(landmarks_c1 .^ 2, 1), 2) );
landmarks_s = landmarks_c1 ./ scales;
landmarks_out = zeros(size(landmarks_s));
angles = zeros(size(landmarks_s, 3), 1);

if isfield(options, 'use_given_rotation_matrix') && options.use_given_rotation_matrix && isfield(options, 'angles')
	angles = options.angles; 
end

for i = 1 : size(landmarks, 3)
    cur_landmark = landmarks_s(:, :, i);
	
	if isfield(options, 'use_given_rotation_matrix') && options.use_given_rotation_matrix
		cur_angle = angles(i);
		rotmat = [cos(cur_angle), -sin(cur_angle); sin(cur_angle), cos(cur_angle)];
		landmark_rot = cur_landmark * rotmat';
	else 
		% perform alignment based on major axis 
		[coeff, score] = pca(cur_landmark);

		cur_angle = -atan2(coeff(2, 1), coeff(1, 1));
		cur_angle = cur_angle + pi / 4;
		rotmat = [cos(cur_angle), -sin(cur_angle); sin(cur_angle), cos(cur_angle)];
		landmark_rot = cur_landmark * rotmat';
		cur_skew = skewness(landmark_rot);
		if sum(cur_skew) < 0
			landmark_rot = landmark_rot * [-1, 0; 0, -1];
			cur_angle = cur_angle + pi;
		end
		cur_angle = (cur_angle > pi) * (-2 * pi) + (cur_angle < -pi) * (2 * pi) + cur_angle;
		angles(i) = cur_angle;
	end

    % 05/25/2018
    if strcmp(component, 'cell') || strcmp(component, 'nuc')
        landmark_shift = shift_landmarks_sequence(landmark_rot);
        landmarks_out(:, :, i) = landmark_shift;
    elseif strcmp(component, 'both')
        landmark_cell_shift = landmark_rot(1 : N_cell_points, :);
        landmark_cell_shift = shift_landmarks_sequence(landmark_cell_shift);

        landmark_nuc_shift = landmark_rot(N_cell_points + 1 : end, :);
        landmark_nuc_shift = shift_landmarks_sequence(landmark_nuc_shift);
        landmarks_out(:, :, i) = [landmark_cell_shift; landmark_nuc_shift];
    end

    if false
        plot(landmark_rot(:, 1), landmark_rot(:, 2), '.')
        xlim([-0.05, 0.05]);
        ylim([-0.05, 0.05]);
    end
end

end


function [landmarks_out] = shift_landmarks_sequence(landmarks_in)
% 05/27/2018 circshift landmarks via alignment with cosine curve.

N_points = size(landmarks_in, 1);
cos_wv = cos(((1 : N_points) / N_points - 1 / 2) * 2 * pi);

cur_xcorr = cconv(cos_wv, landmarks_in(:, 2), N_points);
[~, max_ind] = max(cur_xcorr);
landmarks_out = circshift(landmarks_in, N_points - max_ind + 1);

end

