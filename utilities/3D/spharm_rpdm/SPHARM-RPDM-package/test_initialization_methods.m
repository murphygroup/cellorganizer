function [obj_val_mat] = test_initialization_methods(vertices, faces)
% compute ordered neighbors of all vertices. 

[neighbors, ~, face_ids] = compute_neigbhors(faces, 2);
total_nb_num = sum(cellfun(@numel, neighbors));

% initial parameterization using the one from spharm-mat
% m_x = start_values_function(vertices, faces, neighbors);
[sph_verts] = initParamCALD_x1(vertices, faces);
[obj_val] = goal_func(sph_verts, neighbors, total_nb_num)
[sph_verts_pca] = initParamCALD_x(vertices, faces);
[obj_pca] = goal_func(sph_verts_pca, neighbors, total_nb_num)
[cartesian] = start_values_function(vertices, faces, neighbors);
[obj_orig] = goal_func(cartesian, neighbors, total_nb_num)

%
obj_val_mat = [obj_val, obj_pca, obj_orig];

% compute hausdorff distances
[V, tri] = SpiralSampleSphere(size(vertices, 1));

[hd, lhd, rhd] = hausdorff_distance(V, sph_verts);
[hd_pca, lhd_pca, rhd_pca] = hausdorff_distance(V, sph_verts_pca);

[hd_c, lhd_c, rhd_c] = hausdorff_distance(V, cartesian);

[is_good, hd_r] = spherical_parameterization_quality_check(sph_verts, vertices, faces);
[is_good_1, hd_r1] = spherical_parameterization_quality_check(sph_verts_pca, vertices, faces);
[is_good_2, hd_r2] = spherical_parameterization_quality_check(cartesian, vertices, faces);


end


function [obj] = goal_func(m_x, neighbors, total_nb_num)

nvert = size(m_x, 1);
nbsum_vec = zeros(nvert, 3);
for i = 1 : nvert
    nbsum_vec(i, :) = sum(m_x(neighbors{i}, :));    
end
% obj = sum(m_x .* nbsum_vec, 2);
% obj = 0.5 * sum(cellfun(@numel, neighbors) - obj(:));
obj = 0.5 * (total_nb_num - sum((m_x(:) .* nbsum_vec(:))));

end