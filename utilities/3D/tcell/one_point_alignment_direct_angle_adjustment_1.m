function [opt_angles] = one_point_alignment_direct_angle_adjustment_1(I_seg, I_raw, synapse_center, param)
%
% idea: first get the region close to the center of the synapse by a
% cylinder, then threshold the intensity of the selected region. After that
% for each slice, use weighted least square fitting to decide the line to
% fit the region, the weights are intensities of the pixels. At last, take
% the average of the slope using a Guassian function as weight from center
% slice to two sides. And calculate the angle by z-axis. 
%
%
% Author: Xiongtao Ruan
% Date: Sept. 09, 2016


close all
theta_1 = 0;
theta_2 = 0;
thickness = 15;
shift = -1;
[cy, cx, cz] = ind2sub(size(I_seg), find(I_seg > 0));
coordinates = [cx, cy, cz];

% figure, imshow(reshape_2d(I_raw .* (I_seg >= 0.5)), []);
region_coordinates = generate_region(theta_1, theta_2, thickness, shift, synapse_center, coordinates, param); 

if false
    figure, 
    image_to_show = [];
    image_to_show = [image_to_show; I_seg ./ 2];
    image_to_show(sub2ind(size(I_seg), region_coordinates(:, 2), region_coordinates(:, 1), region_coordinates(:, 3))) = 0.75;    
    image_to_show(34 : 38, 34 : 38, 10 : 26) = 1;    
    I_raw_shown = mat2gray(I_raw .* (I_seg >= 0.5) + (I_seg < 0.5) .* min(I_raw(I_seg >= 0.5))) .* 0.5;
    I_raw_shown(sub2ind(size(I_seg), region_coordinates(:, 2), region_coordinates(:, 1), region_coordinates(:, 3))) = 0.75;
    I_raw_shown(34 : 38, 34 : 38, 10 : 26) = 1;   
    image_to_show = [image_to_show; I_raw_shown];
    imshow(reshape_2d(image_to_show), []);
end

I_seg_region = zeros(size(I_seg));
I_seg_region(sub2ind(size(I_seg), region_coordinates(:, 2), region_coordinates(:, 1), region_coordinates(:, 3))) = 1;

I_raw_region = I_raw .* I_seg_region;
intensity_percent_thrsh = 80;
thresh = prctile(I_raw_region(I_seg_region > 0), intensity_percent_thrsh);

I_seg_reg_thresh = I_raw_region >= thresh;

I_weight = I_raw .* I_seg_reg_thresh;

[~, ~, sz] = size(I_raw);
betas = zeros(sz, 2);
good_inds = zeros(sz, 1);
% threshold of the length in x direction
I_seg_y = squeeze(sum(I_seg_reg_thresh, 1));
x_length = arrayfun(@(z) find(I_seg_y(:, z) > 0, 1, 'last') - find(I_seg_y(:, z) > 0, 1, 'first') + 1, 1 : sz, 'uniformoutput', false);
x_length = cell2mat(x_length(~cellfun(@isempty, x_length)));
x_thrsh_percent = 50;
x_thrsh_length = prctile(x_length, x_thrsh_percent);

for i = 1 : sz
    slice_i = I_weight(:, :, i);
    I_seg_thresh_i = I_seg_reg_thresh(:, :, i);
    inds = find(I_seg_thresh_i > 0);
    
    % set weight as the intensity, because it will counts more for higher
    % intensities.     
    weights = slice_i(inds);
    [y, x] = ind2sub(size(I_seg_thresh_i), inds);

    if isempty(x) || max(x) - min(x) + 1 < x_thrsh_length
        continue;
    end
    W = diag(weights);
    X = [x, ones(size(x, 1), 1)];
    beta_i = (X' * W * X) \ (X' * W * y);
    betas(i, :) = beta_i;
    good_inds(i) = 1;
end

theta_weight = normpdf(1 : sz, floor((sz + 1) / 2), (sz - 1) / 6 / 2);

if false
    figure, imshow(reshape_2d(I_seg_reg_thresh .* I_raw + (I_seg_reg_thresh == 0) .* min(I_raw(I_seg_reg_thresh > 0))), [])
end

beta_mean = theta_weight(good_inds == 1) * betas(good_inds == 1, 1) ./ sum(theta_weight(good_inds == 1));

theta_2 = 0;
theta_2 = atan(beta_mean);

opt_angles = [0, theta_2];

end

function region_coordinates = generate_region(theta_1, theta_2, thickness, cur_shift, synapse_center, coordinates, param)
    % we don't rotation through y-axis
    [M] = euler_angles_to_tranformation_matrix(theta_1, 0, theta_2, 1, param.angle_order);
    v = M * ([0, -1, 0])';
    cur_center = synapse_center + [cur_shift, 0, 0];
    half_thickness = (thickness - 1) / 2;
    inds = abs(bsxfun(@minus, coordinates, cur_center) * v) - half_thickness <= 0;
    region_coordinates = coordinates(inds, :);    
end
