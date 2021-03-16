function [Point]=Trans_Init_1( GFP_raw_seg, GFP_seg, centroid, centroid_apc, left_right_end_points)

L_x = left_right_end_points(1);
L_y = left_right_end_points(2);
R_x = left_right_end_points(3);
R_y = left_right_end_points(4);

line_vector = [L_x - R_x, L_y - R_y];
line_vector = line_vector / norm(line_vector, 2);
line_norm_vector = [line_vector(2), -line_vector(1)];
B = bwboundaries(GFP_seg);
bw_coordinates = B{1};

% interpolate the boundary coordinates

x_boundary = bw_coordinates(:, 2);
y_boundary = bw_coordinates(:, 1);
% smooth the boundary points
x_boundary = smooth(x_boundary);
y_boundary = smooth(y_boundary);

x_boundary = [x_boundary; x_boundary(1)];
y_boundary = [y_boundary; y_boundary(1)];

inds = 1 : numel(x_boundary);
q_inds = 1 : 0.1 : numel(x_boundary);

qx_boundary = interp1(inds, x_boundary, q_inds);
qy_boundary = interp1(inds, y_boundary, q_inds);
boundary_vec = [qx_boundary', qy_boundary'];

% find the nearest points to the line containing the centroid of apc cell
% 
distances = (boundary_vec - centroid_apc) * line_norm_vector';

% Centroid_frame_0 = (left_right_end_points(1:2) + left_right_end_points(3:4)) / 2;
% Centroid_far = (Centroid_frame_0 - centroid) * 2 + centroid;
% distances = (boundary_vec - Centroid_far) * line_norm_vector';

if ~(all(distances >= 0) || all(distances <= 0))
    % in this case, it means the line containing centroid of apc cell has
    % intersection point with the cell boundary.
    Center_frame_0 = (left_right_end_points(1:2) + left_right_end_points(3:4)) / 2;
    Centroid_far = (Center_frame_0 - centroid) * 3 + centroid;
    distances = (boundary_vec - Centroid_far) * line_norm_vector';
    % error('special case, not implement yet!')   
end

if all(distances >= 0)
    line_norm_vector = line_norm_vector / norm(line_norm_vector, 2);
elseif all(distances <= 0)
    line_norm_vector = -line_norm_vector / norm(line_norm_vector, 2);
end

[~, nearest_ind] = min(abs(distances));
nearest_coordinates = boundary_vec(nearest_ind, :);

% thresh_coord = nearest_coordinates * 0.8 + centroid * 0.2;
thresh_coord = nearest_coordinates * 2 / 3 + centroid * 1 / 3;
hard_thresh_depth_pixel = 5;

thresh_depth = min(hard_thresh_depth_pixel, (thresh_coord - nearest_coordinates) * line_norm_vector');

depths = linspace(1, thresh_depth, 100);

dist_1 = (centroid - [L_x + R_x, L_y + R_y] / 2) * line_norm_vector';
dist_2 = (centroid - nearest_coordinates) * line_norm_vector';
if dist_1 / dist_2 > 0.7
    Point = [];
    return;
end

best_L_i = [];
best_R_i = [];
best_intensity = 0;
best_lateral_dist = 0;
for i = 1 : numel(depths)
    depth_i = depths(i);
    [L_i, R_i, intensity_i] = extract_proporty_function(line_norm_vector, nearest_coordinates, depth_i, boundary_vec, GFP_raw_seg);
    centroid_lateral_dist = line_vector * ((L_i + R_i) / 2 - [L_x + R_x, L_y + R_y] / 2)';
    if intensity_i > best_intensity && pdist2(L_i, R_i) > 5 && pdist2(L_i, R_i) > 0.75 * pdist2([L_x, L_y], [R_x, R_y]) && abs(centroid_lateral_dist) < 0.2 * pdist2([L_x, L_y], [R_x, R_y])
        best_L_i = L_i;
        best_R_i = R_i;
        best_intensity = intensity_i;
        best_lateral_dist = centroid_lateral_dist;
    end
end
if false && ~isempty(best_L_i) && ~isempty(best_R_i)
    figure(1), imshow(GFP_raw_seg);
    hold on, plot(best_L_i(1), best_L_i(2), 'r*', best_R_i(1), best_R_i(2), 'r*');
    hold on, plot(L_x, L_y, 'ro', R_x, R_y, 'ro');
    hold off
    figure(2), plot(boundary_vec(:, 1), boundary_vec(:, 2));
    hold on, plot(best_L_i(1), best_L_i(2), 'r*', best_R_i(1), best_R_i(2), 'r*');
    hold on, plot(L_x, L_y, 'ro', R_x, R_y, 'ro');
    hold off
end

Point = [best_L_i, best_R_i];
% best_lateral_dist

end


function [L_i, R_i, intensity_i] = extract_proporty_function(line_norm_vector, nearest_coordinates, depth, boundary_vec, GFP_raw_seg)
    x0 = depth * line_norm_vector + nearest_coordinates;
    bound_dists = abs((boundary_vec - x0) * line_norm_vector');
    
    [~, nearest_ind] = min(bound_dists);
    L_i = boundary_vec(nearest_ind, :);
    
    [~, inds] = sort(bound_dists);
    chosen_coords = boundary_vec(inds(2:10), :);
    hard_thesh = 0.5;
    
    L_dists = pdist2(L_i, chosen_coords);
    [~, farthest_ind] = max(L_dists);
    R_i = chosen_coords(farthest_ind, :);
    for i = 1 : size(chosen_coords, 1)
        cur_coord = chosen_coords(i, :);
        cur_L_dist = pdist2(L_i, cur_coord);
        if cur_L_dist > hard_thesh || cur_L_dist > pdist2(x0, L_i)
            R_i = cur_coord;
            break
        end         
    end
    
    % get the mean intensity across the line
    
    Center = (L_i + R_i) / 2;
    half_dist_thresh = 2;
    if pdist2(Center, L_i) < half_dist_thresh
        half_dist_thresh =  pdist2(Center, L_i);
    end
    line_vector = [line_norm_vector(2), -line_norm_vector(1)];
    
    points_in_line = Center + (-half_dist_thresh : 0.01 : half_dist_thresh)' * line_vector;
    
    [X, Y] = meshgrid(1 : size(GFP_raw_seg, 2), 1 : size(GFP_raw_seg, 1));
    intensities = interp2(X, Y, GFP_raw_seg, points_in_line(:, 1), points_in_line(:, 2));
    intensity_i = mean(intensities);

end


