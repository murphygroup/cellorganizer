function [is_shifting, bad_timepoint_inds] = synapse_shifting_detection_function(synapses_mat, param)
% detect whether synapse shift between different cells. 

if nargin < 1
    load('test_tcell.mat');
end

shift_len_thrsh = 15;
shift_mov_len_thrsh = 12.5;
cos_shift_angle_thrsh = -150;
shift_num_thrsh = 2;

synapse_coords = synapses_mat(:, 1 : 4);
synapse_centers = (synapse_coords(:, 1 : 2) + synapse_coords(:, 3 : 4)) / 2;

synapse_shift_vec = diff(synapse_centers);

synapse_shift_len = sqrt(sum(synapse_shift_vec .^ 2, 2));

mov_diff = synapse_centers - movmedian(synapse_centers, 3);

mov_diff_len = sqrt(sum(mov_diff .^ 2, 2));


cos_shift_angle = sum(synapse_shift_vec(1 : end - 1, :) .* synapse_shift_vec(2 : end, :), 2);
cos_shift_angle = [1; cos_shift_angle];

judge_mat = [synapse_shift_len > shift_len_thrsh, mov_diff_len(2 : end) > shift_mov_len_thrsh, cos_shift_angle < cos_shift_angle_thrsh];

inds = find(sum(judge_mat, 2) >= 2);

% remove those suspected time point and see if there is still dramatic
% change to decide if there is shift. 

bad_timepoint_mat = false(numel(inds), 1);

for i = 1 : numel(inds)
    synapse_centers_i = synapse_centers;
    cur_ind = inds(i);
    synapse_centers_i(cur_ind, :) = [];
    
    synapse_shift_vec_i = diff(synapse_centers_i);
    synapse_shift_len_i = sqrt(sum(synapse_shift_vec_i .^ 2, 2));
    
    mov_diff_i = synapse_centers_i - movmedian(synapse_centers_i, 3);
    mov_diff_len_i = sqrt(sum(mov_diff_i .^ 2, 2));
    
    % cos_shift_angle_i = sum(synapse_shift_vec_i(1 : end - 1, :) .* synapse_shift_vec_i(2 : end, :), 2);
    % cos_shift_angle_i = [1; cos_shift_angle_i];

    judge_mat_i = [synapse_shift_len_i > shift_len_thrsh, mov_diff_len_i(2 : end) > shift_mov_len_thrsh];
    
    if sum(judge_mat_i(cur_ind, :)) <= 0
        bad_timepoint_mat(i) = true;
    end
end

bad_timepoint_inds = inds(bad_timepoint_mat);

is_shifting = false;
if sum(bad_timepoint_mat) > shift_num_thrsh
    is_shifting = true;
end



end



