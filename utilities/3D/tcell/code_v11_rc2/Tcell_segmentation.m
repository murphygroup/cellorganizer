function [loc, I_5_final] = Tcell_segmentation(I_in)
% use local maximum and local minimum based segmentation
%
% Author: Xiongtao Ruan
%


close all;
debug = false;

if nargin < 1
    load('GFP_v2.mat');
    I_in = GFP{15};
    debug = ~true;
    % I_in = GFP{13};
end

I_1 = imgaussfilt(mat2gray(I_in), 0.8);

I_2 = I_1;

I_2_1 = imtophat(I_2,strel('disk',15));

I_2_2= imadjust(I_2_1);

I_2 = I_2_2;

thrsh = graythresh(I_2 );

background = I_2(I_2 <= 0.5 * thrsh);

BW_min = imregionalmin(I_2);

%I_3 = imimposemin(I_2, BW_min);

%L = watershed(I_3);

I_4 = I_2 > max(background);
I_4_1 = imfill(I_4, 'hole');
I_4_2 = imopen(I_4_1, strel('disk', 5));
I_4_3 = imerode(I_4_2, strel('disk', 3));
% I_4_4 = imerode(I_4_2, strel('line', 7, 0));
% I_4_5 = imerode(I_4_4, strel('line', 7, 90));

BW_min_1 = I_4_3 .* BW_min;

I_5 = imimposemin(I_2, BW_min_1);

L_1 = watershed(I_5);
L_1 = double(L_1);

if debug
    figure, imshow(L_1, []);
end

% BW_max = imregionalmax(I_2);
% if debug
%     figure, imshow(BW_max);
% end
% 
% BW_max_1 = I_4_2 .* BW_max;
% 
% I_6 = imimposemin(I_2, BW_max_1);
% 
% L_2 = watershed(I_6);
% if debug
%     figure, imshow(L_2, []);
% end

background_label = mode(L_1(:) .* (1 - I_4_2(:)));

L_3 = L_1;

L_3(L_1 == background_label) = 0;

se = strel('disk', 1);

L_3_1 = imclose(L_3, se);

% se = strel('disk', 1);
% 
% L_3_2 = imerode(L_3_1, se);


L_3_labels = bwlabel(L_3_1 > 0);
% all_stats = regionprops(L_3_1 > 0, 'Area', 'Eccentricity', 'Solidity', 'MajorAxisLength', 'MinorAxisLength');
all_stats = regionprops(L_3_1 > 0, 'Area', 'MajorAxisLength');
All_area = [all_stats.Area];
All_major_axis = [all_stats.MajorAxisLength];

prctile_low = prctile(All_area, 10);
prctile_high = prctile(All_area, 90);
mean_area = mean(All_area(All_area > prctile_low & All_area < prctile_high));
std_area = std(All_area(All_area > prctile_low & All_area < prctile_high));

major_axis_low = prctile(All_major_axis, 10);
major_axis_high = prctile(All_major_axis, 90);
mean_major_axis = mean(All_major_axis(All_major_axis > major_axis_low & All_major_axis < major_axis_high));
std_major_axis = std(All_major_axis(All_major_axis > major_axis_low & All_major_axis < major_axis_high));

labels = unique(L_3_labels(:));
labels(labels == 0) = [];

L_3_final = L_3;

for i = 1 : numel(labels)
    cur_label = labels(i);
    
    L_3_i = L_3_labels == cur_label;
    
    L_3_parts = L_3 .* L_3_i;
    if debug
        figure(1), imshow(L_3_parts, [])
    end
    
    BW_min_i = BW_min_1 .* L_3_i;
    BW_neighbor = bwdist(L_3_i) < 3;
%     BW_max_neighbor = BW_max_1 .* BW_neighbor;
  
    if mean(I_2(BW_neighbor == 1) ) > mean(background)  + 2 * std(background)
        is_cell = true;
    else
        %%%L_3_final(L_3_i == 0);%%%jing dan yang gai
        L_3_final(L_3_i == 1) = 0;
        disp('There is no valid cell here, skip it!');
        continue;
    end
    
%     local_max_inds = find(BW_max_neighbor == 1);
%     [lmax_y, lmax_x] = ind2sub(size(BW_max_neighbor), local_max_inds);
%     lmax_coords = [lmax_x, lmax_y];
    
    local_min_inds = find(BW_min_i == 1);
    [lmin_y, lmin_x] = ind2sub(size(BW_min_i), local_min_inds);
    lmin_coords = [lmin_x, lmin_y];

    % size(lmin_coords, 1)
    % size(lmax_coords, 1)
    
    if size(lmin_coords, 1) <= 5
        
        % if there is a single local minimum, then check if there are other
        % local minimum surround, if it is, then just remove it.
        
        if size(lmin_coords, 1) == 1
            if sum(L_3_i(:)) < 25
                L_3_i_bw_surround = bwdist(L_3_i);
                surround_img = ((L_3_i_bw_surround > 0) & (L_3_i_bw_surround < 10)) .* L_3_labels;
                if any(surround_img(:))
                    L_3_final(L_3_i == 1) = 0;
                    continue;
                end
                
            end
            
        end
        
%         lmin_max_dists = pdist2(lmin_coords, lmax_coords);
%         
%         ratios = max(lmin_max_dists, [], 2) ./ min(lmin_max_dists, [], 2)      
%         
%         [~, best_ind] = min(max(ratios, [], 2));
%         chosen_seed = lmin_coords(best_ind, :);

        % cur_stats = regionprops(L_3_i > 0, 'Area', 'Eccentricity', 'Solidity', 'MajorAxisLength', 'MinorAxisLength');
        cur_stats = regionprops(L_3_i > 0, 'Area', 'Eccentricity','MajorAxisLength');
                
        is_multiple_cell = false;
        if cur_stats.Area > mean_area + 3 * std_area && cur_stats.MajorAxisLength > mean_major_axis + 4 * std_major_axis
            is_multiple_cell = true;
        end
        if cur_stats.Area > mean_area + 2 * std_area && cur_stats.MajorAxisLength > mean_major_axis + 5 * std_major_axis
            is_multiple_cell = true;
        end
        if cur_stats.Area > mean_area + 1 * std_area && cur_stats.MajorAxisLength > mean_major_axis + 3 * std_major_axis && cur_stats.Eccentricity > 0.95
            is_multiple_cell = true;
        end

        if   is_multiple_cell
            
            disp('There are more than one cell, and most probably two cells, use the two farest local minimal')
            lmin_dists = pdist(lmin_coords);
            dist_mat = squareform(lmin_dists);
            [~, max_ind] = max(dist_mat(:));
            [ind_1, ind_2] = ind2sub(size(dist_mat), max_ind);
            part_label_1 = L_3_parts(lmin_coords(ind_1, 2), lmin_coords(ind_1, 1));
            part_label_2 = L_3_parts(lmin_coords(ind_2, 2), lmin_coords(ind_2, 1));
            L_3_final(L_3_i == 1) = 0;
            L_3_parts_1 = L_3_parts == part_label_1 | L_3_parts == part_label_2;
            L_3_parts_1 = imopen(L_3_parts_1, strel('disk', 3));
            L_3_final(L_3_parts_1) = ceil((part_label_1 + part_label_2) / 2);
        end
    else
        
        if any(pdist(lmin_coords) < 10)
            is_cell = true;
        end
        
%         lmin_max_dists = pdist2(lmin_coords, lmax_coords);
%         
%         ratios = max(lmin_max_dists, [], 2) ./ min(lmin_max_dists, [], 2)    
        
        % cur_stats = regionprops(L_3_i > 0, 'Area', 'Eccentricity', 'Solidity', 'MajorAxisLength', 'MinorAxisLength');
        cur_stats = regionprops(L_3_i > 0, 'Area', 'Eccentricity', 'MajorAxisLength');
        
        is_multiple_cell = false;
        if cur_stats.Area > mean_area + 3 * std_area
            is_multiple_cell = true;
        end
        if cur_stats.Area > mean_area + 2 * std_area && cur_stats.MajorAxisLength > mean_major_axis + 3 * std_major_axis
            is_multiple_cell = true;
        end
        if cur_stats.Area > mean_area + 1.5 * std_area && cur_stats.MajorAxisLength > mean_major_axis + 3 * std_major_axis && cur_stats.Eccentricity > 0.92
            is_multiple_cell = true;
        end

        if   is_multiple_cell
            % there are multiple cells
            disp('probably two cells')
            
            % separate them by try to remove any part and see if it can be
            % broken to two parts, which has similar area
            
            cur_unique_inds = unique(L_3_parts);
            cur_unique_inds(cur_unique_inds == 0) = [];
            sep_areas = [];
            for j = 1 : numel(cur_unique_inds)
                I_part_j = L_3_parts == cur_unique_inds(j);
                
                I_part_j_1 = imdilate(I_part_j, strel('disk', 1));
                L_3_sep_j = L_3_i - I_part_j_1;
                L_3_sep_j(L_3_sep_j < 0) = 0;
                CC = bwconncomp(L_3_sep_j);
                if CC.NumObjects == 2
                    cc_area = cellfun(@numel, CC.PixelIdxList);
                    sep_areas = [sep_areas; [j, cc_area]];          
                elseif CC.NumObjects > 2
                    cc_area = cellfun(@numel, CC.PixelIdxList);
                    cc_area = sort(cc_area);
                    sep_areas = [sep_areas; [j, cc_area(end-1:end)]];                            
                end
            end
            
            if isempty(sep_areas)
                disp('There are more than one cell, and most probably two cells, use the two farest local minimal')
                lmin_dists = pdist(lmin_coords);
                dist_mat = squareform(lmin_dists);
                [~, max_ind] = max(dist_mat(:));
                [ind_1, ind_2] = ind2sub(size(dist_mat), max_ind);
                part_label_1 = L_3_parts(lmin_coords(ind_1, 2), lmin_coords(ind_1, 1));
                part_label_2 = L_3_parts(lmin_coords(ind_2, 2), lmin_coords(ind_2, 1));
                L_3_final(L_3_i == 1) = 0;
                L_3_parts_1 = L_3_parts == part_label_1 | L_3_parts == part_label_2;
                L_3_parts_1 = imopen(L_3_parts_1, strel('disk', 3));
                L_3_final(L_3_parts_1) = ceil((part_label_1 + part_label_2) / 2);                        
            else         
                [~, ind] = min(abs(diff(sep_areas(:, 2:3), [], 2)));
                L_3_final(L_3_final == cur_unique_inds(sep_areas(ind, 1))) = 0; 
            end
        end
    end 
    
    if debug
        figure(2), imshow(L_3_final, [])
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_3_final = imclose(L_3_final > 0, strel('disk', 1));
% I_3_final_bound = I_3_final;
% I_3_final_bound(1 : end, 1) = 1;
% I_3_final_bound(1 : end, end) = 1;
% I_3_final_bound(1, 1 : end) = 1;
% I_3_final_bound(end, 1: end) = 1;
% 
% I_3_seg = imimposemin(I_2, I_3_final_bound);
% 
% L_3 = watershed(I_3_seg);
% L_3 = double(L_3);
% 
% figure, imshow(L_3);

% in some case some very bright cell has no local minimum, in this case, we
% judge, if there are high intensity object not labeled, then add it to the
% labeling. 

I_2_1 = imtophat(I_1,strel('disk',15));

I_2_2= imadjust(I_2_1);
background_1 = I_2_2(I_2_2 <= 0.4 * thrsh);
% BW_min = imregionalmin(I_2_2);

I_5 = I_2_2 > max(background_1);
I_5_1 = imfill(I_5, 'hole');
I_5_2 = imopen(I_5_1, strel('disk', 5));
I_5_3 = imerode(I_5_2, strel('disk', 3));

I_5_4 = imclearborder(I_5_3);
CC_1 = bwconncomp(I_5_4);

for i = 1 : CC_1.NumObjects
    I_obj = zeros(size(I_3_final));
    cur_pixel_list = CC_1.PixelIdxList{i};
    I_obj(cur_pixel_list) = 1;
    if ~I_3_final(cur_pixel_list) > 0        %
        BW_neighbor = bwdist(I_obj) < 2;  
        if mean(I_2(BW_neighbor == 1) ) > mean(background_1)  + 1.5 * std(background_1)
            is_cell = true;
        else
            disp('There is no valid cell here, skip it!');
            continue;
        end
        % figure, imshow(I_obj, []);
        I_3_final(cur_pixel_list) = 1;      
    end
end

stats = regionprops(I_3_final);

loc = cat(1, stats.Centroid);

% make a combination of the old method and this method, in case of missing
% some cells with low intensity. 

loc_old = MarkerControlled_Watershed2(I_in);
loc_old = cat(1, loc_old{:});

all_dists = pdist2(loc_old, loc);
min_dists = min(all_dists, [], 2);
min_dist_low = prctile(min_dists, 10);
min_dist_high = prctile(min_dists, 90);
min_dists_truncate = min_dists(min_dists > min_dist_low & min_dists < min_dist_high);
min_dist_mean = mean(min_dists_truncate);
min_dist_std = std(min_dists_truncate);

miss_inds = find(min_dists > min_dist_mean + 4 * min_dist_std);

for i = 1 : numel(miss_inds)
    loc_i = round(loc_old(miss_inds(i), :));
%     dist_i = pdist2(loc_i, loc);
%     figure, imshow(I_2, []);
%     hold on, 
%     plot(loc_i(:, 1), loc_i(:, 2), 'o')
%     hold off 
    
    I_obj = zeros(size(I_3_final));
    I_obj(loc_i(2), loc_i(1)) = 1;
    BW_neighbor = bwdist(I_obj) < 5;  
    if mean(I_2(BW_neighbor == 1) ) > mean(background_1)  + 1.5 * std(background_1)
        is_cell = true;
    else
        disp('There is no valid cell here, skip it!');
        continue;
    end
    I_3_final(I_obj == 1) = 1;    
end

stats = regionprops(I_3_final);
loc = cat(1, stats.Centroid);

if ~false && debug
    figure, imshow(I_2, []);
    hold on, 
    plot(loc(:, 1), loc(:, 2), 'ro')
    hold on, 
    plot(loc_old(:, 1), loc_old(:, 2), 'b*')
    hold off
    figure, imshow(L_3_final, [])
    hold on, 
    plot(loc(:, 1), loc(:, 2), 'ro')
    hold off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try to segment the cells after getting the makers. using intensity based
% method, and find corners to separate touching cells. 

% I_2_1 = imerode(I_2, strel('disk', 3));
% 
% mask = zeros(size(I_2));
% mask(sub2ind(size(I_2), round(loc(:, 2)), round(loc(:, 1)))) = 1;
% mask = imdilate(mask, strel('disk', 3));
% bw = activecontour(I_2,mask);

% I_6 = imimposemin(imcomplement(I_2_1), mask);
% L_2 = watershed(I_6);

% I_2_1 = imtophat(I_1,strel('disk',15));
thrsh = graythresh(I_2);

background_1 = I_2_2(I_2_2 <= 0.5 * thrsh);
I_5 = I_2_2 > max(background_1);
I_5_1 = imfill(I_5, 'hole');
I_5_final = I_5_1;

I_5_label = bwlabel(I_5_1);
unique_labels = unique(I_5_label(:));
unique_labels(unique_labels == 0) = [];

I_obj = zeros(size(I_5_1));
loc_1 = round(loc);
I_obj(sub2ind(size(I_obj), loc_1(:, 2), loc_1(:, 1))) = 1;

for i = 1 : numel(unique_labels)
    cur_label = unique_labels(i);
    I_cur_cell = I_5_label == cur_label;
    num_cell = sum(sum((I_cur_cell .* I_obj)));
    if num_cell < 2
        continue;
    end
    loc_inds = find(I_cur_cell .* I_obj > 0);
    [loc_y, loc_x] = ind2sub(size(I_obj), loc_inds);
    loc_points = [loc_x, loc_y];
    
    % merge very close centers, which should be for the same cell. 
    if size(loc_points, 1) >= 2
        loc_points_1 = loc_points;
        loc_dists = pdist(loc_points);
        if any(loc_dists < 15)
            loc_dists_1 = squareform(loc_dists);
            loc_dists_2 = tril(loc_dists_1, -1);
            loc_dists_2(loc_dists_2 == 0) = 10000;
            loc_dist_inds = find(loc_dists_2 < 15);
            for j = 1 : numel(loc_dist_inds)
                [loc_ind_1, loc_ind_2] = ind2sub(size(loc_dists_1), loc_dist_inds(j));
                loc_points_1([loc_ind_1, loc_ind_2], :) = [];
                loc_points_1 = round([loc_points_1; (loc_points(loc_ind_1, :) + loc_points(loc_ind_2, :)) / 2]);
                break;
            end
            loc_points = loc_points_1;
            num_cell = num_cell - 1;
        end
        if num_cell < 2
            continue;
        end
    end
    corner = detectMinEigenFeatures(I_cur_cell);
    if debug
        figure, imshow(I_cur_cell, []);
        hold on, plot(corner.selectStrongest(num_cell * 4));
    end
%     I_cur_loc = 1;   

    corner_1 = corner.selectStrongest(num_cell * 4);
    corner_points = double(corner_1.Location);
    
    boundary = extract_cell_boundary(I_cur_cell);
    
    % match corner points to the boundary
    corner_boundary_dists = pdist2(corner_points, boundary);
    [~, inds] = min(corner_boundary_dists, [], 2);
    
    mapped_corner_points = boundary(inds, :);
    
    eu_dists = pdist2(mapped_corner_points, mapped_corner_points);
    
    all_pair_dists = sqrt(sum((boundary - [boundary(end, :); boundary(1:end-1, :)]) .^ 2, 2));
    curcumstance = sum(all_pair_dists);
    
    all_cum_dists = cumsum(all_pair_dists);
%     plot(loc(:, 1), loc(:, 2), 'ro')

    dists_to_zero = all_cum_dists(inds);
    pair_curve_dists = dists_to_zero - dists_to_zero';
    pair_curve_dists_1 = abs(pair_curve_dists);
    pair_curve_dists_2 = curcumstance - abs(pair_curve_dists);
    pair_curve_dists_f = min(pair_curve_dists_1, pair_curve_dists_2);
    
    dist_ratios = eu_dists ./ pair_curve_dists_f;
    
    dist_ratios_l = tril(dist_ratios, -1);
    dist_ratios_l(dist_ratios_l == 0) = 1;
    dist_ratios_array = dist_ratios_l(:);
    [sorted_ratios, sort_inds] = sort(dist_ratios_array);
    
    % addpath('./code_cell_pair_v0/bresenham_2d');
    all_lines = {};
    loc_groups = {};
    for j = 1 : num_cell * 6
%         cur_ratio = sorted_ratios(j, :);
        cur_ind = sort_inds(j);
        [cur_ind_1, cur_ind_2] = ind2sub(size(dist_ratios), cur_ind);
        point_1 = mapped_corner_points(cur_ind_1, :);
        point_2 = mapped_corner_points(cur_ind_2, :);
        point_dist = sqrt(sum((point_1 - point_2) .^ 2)); 
        if point_dist > 20 || pair_curve_dists_f(cur_ind_1, cur_ind_2) < 25
            continue;
        end
        [x, y]=bresenham(point_1(1),point_1(2),point_2(1),point_2(2));
        I_line = zeros(size(I_5_1));
        I_line(sub2ind(size(I_line), y, x)) = 1;
        I_line_1 = imdilate(I_line, strel('disk', 1));
        I_rm_line = I_cur_cell;
        I_rm_line(sub2ind(size(I_rm_line), y, x)) = 0;
        if debug
            figure, imshow(I_rm_line - I_line_1, []);        
            hold on, plot(loc_points(:, 1), loc_points(:, 2), 'o');
            hold on, plot(point_1(:, 1), point_1(:, 2), '*', point_2(:, 1), point_2(:, 2), '*');
        end
        if num_cell == 2
            [loc_groups] = group_points_based_on_image(I_cur_cell, point_1, point_2, loc_points);
            if numel(loc_groups) == 2
                I_5_final(I_line_1 == 1) = 0;
                break;       
            end
        end
        
        % for multiple cells, use a greedy method to separate cells one by
        % one
        if num_cell > 2
            % first try to increase threshold to separate the cells            
            if isempty(all_lines)
                all_lines{end+1} = [point_1, point_2];
                % [loc_groups] = group_points_based_on_line(point_1, point_2, loc_points);  
                [loc_groups] = group_points_based_on_image(I_cur_cell, point_1, point_2, loc_points);
            else 
                group_nums = cellfun(@(x) size(x, 1), loc_groups);
                if any(group_nums > 1)
                    group_inds = find(group_nums > 1);
                    
                    for k = 1 : numel(group_inds)
                        cur_group_ind = group_inds(k);
                        cur_points = loc_groups{cur_group_ind};
                        % [cur_groups_k] = group_points_based_on_line(point_1, point_2, cur_points);    
                        [cur_groups_k] = group_points_based_on_image(I_cur_cell, point_1, point_2, cur_points);    
                        if numel(cur_groups_k) >= 2
                            loc_groups(cur_group_ind) = [];
                            % loc_groups = [loc_groups, cur_groups_k];
                            loc_groups = [loc_groups; cur_groups_k];
                            all_lines{end+1} = [point_1, point_2];   
                        end     
                    end
                else
                    break;
                end
            end    
        end
        
    end
    
    % remove lines
    for k = 1 : numel(all_lines)
        cur_points = all_lines{k};
        point_1 = cur_points(1:2);
        point_2 = cur_points(3:4);
        % point_dist = sqrt(sum((point_1 - point_2) .^ 2)); 
        [x, y]=bresenham(point_1(1),point_1(2),point_2(1),point_2(2));
        I_line = zeros(size(I_5_1));
        I_line(sub2ind(size(I_line), y, x)) = 1;
        I_line_1 = imdilate(I_line, strel('disk', 1));
        I_rm_line = I_cur_cell;
        I_rm_line(sub2ind(size(I_rm_line), y, x)) = 0;
        if debug
            figure, imshow(I_rm_line - I_line_1, []);
        end
        I_5_final(I_line_1 == 1) = 0;
    end
end


% future process cells so that get cells with low intensity. 
% to be done.  



% remove small regions that have no loc associated with.
I_5_label = bwlabel(I_5_final);
unique_loc_labels = unique(I_5_label(I_obj > 0));
unique_loc_labels(unique_loc_labels == 0) = [];
I_5_final = ismember(I_5_label, unique_loc_labels);


% remove small regions with area threshold
CC = bwconncomp(I_5_final, 8);
areas = cellfun(@numel, CC.PixelIdxList);
area_low = prctile(areas, 10);
area_high = prctile(areas, 90);
mean_area = mean(areas(areas > area_low & areas < area_high));
std_area = std(areas(areas > area_low & areas < area_high));

area_hard_low = 200;
area_hard_high = 1000;

area_to_rm = areas > min(area_hard_high, mean_area + 3 * std_area) | areas < max(area_hard_low, mean_area - 3 * std_area);
pixels_to_rm = CC.PixelIdxList(area_to_rm);
pixels_to_rm = cat(1, pixels_to_rm{:});
I_5_final(pixels_to_rm) = false;

stats = regionprops(I_5_final);

loc = cat(1, stats.Centroid);

if ~false && debug
    figure, imshow(I_5_final, []);
    hold on, plot(loc(:, 1), loc(:, 2), 'ro')
    hold off
end


end


function [boundary] = extract_cell_boundary(I_in)
% extract boundary points and convert to coordinates. 

boundaries = bwboundaries(I_in);

boundary = flip(boundaries{1}, 2);

end


function [point_groups] = group_points_based_on_line(line_point_1, line_point_2, q_points)

point_signs = sign((q_points - line_point_1) * [line_point_1(2) - line_point_1(2); line_point_2(1) - line_point_1(1)]);

point_groups = {};
if any(point_signs > 0)
    point_groups{end+1} = q_points(point_signs > 0, :);
end
if any(point_signs < 0)
    point_groups{end+1} = q_points(point_signs < 0, :);
end

if any(point_signs == 0)
    point_groups{end+1} = q_points(point_signs == 0, :);
end 

end


function [point_groups] = group_points_based_on_image(I_cell, point_1, point_2, q_points)


[x, y]=bresenham(point_1(1),point_1(2),point_2(1),point_2(2));
I_line = zeros(size(I_cell));
I_line(sub2ind(size(I_line), y, x)) = 1;
I_line_1 = imdilate(I_line, strel('disk', 1));
I_rm_line = I_cell;
I_rm_line = I_rm_line - I_line_1;

I_label = bwlabel(I_rm_line > 0);

q_labels = I_label(sub2ind(size(I_cell), q_points(:, 2), q_points(:, 1)));

unique_labels = unique(q_labels);

point_groups = cell(numel(unique_labels), 1);
for i = 1 : numel(unique_labels)
    point_groups{i} = q_points(q_labels == unique_labels(i), :);
end 

end
