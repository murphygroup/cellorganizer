function [landmarks_out, scales, centers] = convert_to_preshape(landmarks, options)
% convert the landmarks to preshape and the same time save the scale factor and centers. 

% landmarks_c = squeeze(landmarks(:, 1, :) + 1i * landmarks(:, 2, :));

if nargin < 2
    options = struct();
end

options = process_options_structure(struct('use_given_center', false), options);

landmarks_c = squeeze(landmarks(:, 1, :) + 1i * landmarks(:, 2, :));

if isfield(options, 'use_given_center') && options.use_given_center
    centers = options.centers;
else
	centers = mean(landmarks_c, 1);
end
landmarks_c1 = landmarks_c - centers;

scales = sqrt(sum(abs(landmarks_c1) .^ 2));
landmarks_out = landmarks_c1 ./ scales;

end

