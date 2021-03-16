function [hd, lhd, rhd] = hausdorff_distance(X_s, X_t)
% compute hausdorf distance

N_1 = size(X_s, 1);
N_2 = size(X_t, 1);
dist_elem_number_thresh = 5e8;

if N_1 * N_2 > dist_elem_number_thresh
    N_per_ptn = floor(dist_elem_number_thresh / N_1);
    N_ptn = ceil(N_2 / N_per_ptn);
    min_dist_row_mat = 1000 * ones(1, N_1);
    min_dist_col_mat = zeros(N_2, 1);
    for i = 1 : N_ptn
        inds = (i - 1) * N_per_ptn + 1 : min(N_2, i * N_per_ptn);
        X_t_i = X_t(inds, :);
        dist_mat_i = pdist2(X_t_i, X_s);
        cur_min_dist_row_mat = min(dist_mat_i);
        min_dist_row_mat = min([min_dist_row_mat; cur_min_dist_row_mat]);
        min_dist_col_mat(inds) = min(dist_mat_i, [], 2);
    end
    lhd = max(min_dist_col_mat);
    rhd = max(min_dist_row_mat); 
else
    D_mat = pdist2(X_t, X_s);
    lhd = max(min(D_mat, [], 1));
    rhd = max(min(D_mat, [], 2));
end

hd = max(lhd, rhd);

end