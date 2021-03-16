function [is_touching] = touching_tcell_filtering_function(Tcell_segment_cell, Tcell_segment_loc_cell, synapses_mat, param)
% decide if there is neighboring tcell that touches the target one, with a
% distance threshold. 

if nargin < 1
    load('test_tcell.mat');
end

is_touching = false;

for i = 1 : size(synapses_mat, 1)
    cur_frame = synapses_mat(i, end - 1);
	cur_tcell_img = Tcell_segment_cell{cur_frame};
	[cells_dist_mat, L] = calculate_cell_pairwise_distance(cur_tcell_img);
	cur_synapse = synapses_mat(i, :);
	synapse_center = (cur_synapse(1 : 2) + cur_synapse(3 : 4)) / 2;
    I_1 = false(size(cur_tcell_img));
    I_1(round(synapse_center(2)), round(synapse_center(1))) = true;
    se = strel('disk', 3);
    I_1 = imdilate(I_1, se);
    L_labels = L(I_1);
    if ~any(L_labels)
        continue;
    end
    cur_label = mode(L_labels(L_labels ~= 0));
    



end
	




end


function [cells_dist_mat, L] = calculate_cell_pairwise_distance(cur_tcell_img)
% compute pairwise distance between cells

CC = bwconncomp(cur_tcell_img);
L = labelmatrix(CC);

cells_dist_mat = zeros(CC.NumObjects);

cells_subs = cellfun(@(x) custom_ind2sub(size(cur_tcell_img), x), CC.PixelIdxList, 'UniformOutput', false);

for i = 1 : CC.NumObjects
    for j = i + 1 : CC.NumObjects
        cells_dist_mat(i, j) = nearest_distance_two_sets(cells_subs{i}, cells_subs{j});
        cells_dist_mat(j, i) = cells_dist_mat(i, j);
    end
end

end


function [subs] = custom_ind2sub(mat_size, inds)

[y, x] = ind2sub(mat_size, inds);
subs = [x, y];

end


function [hd] = nearest_distance_two_sets(X, Y)

[~, d] = knnsearch(X, Y);
hd = min(d);
    
end
