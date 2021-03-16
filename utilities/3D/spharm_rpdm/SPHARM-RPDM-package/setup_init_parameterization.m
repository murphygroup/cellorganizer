function [m_x] = setup_init_parameterization(vertices, faces, options)
% set up initial parameterization 

% initial parameterization using the one from spharm-mat
% m_x = start_values_function(vertices, faces, neighbors);
% 04/08/2018 for the best one, use hausdorff distances to decide. 


if nargin < 3
    options.initialization_method = 'default';
end

switch options.initialization_method
    case {'default', 'graph_diameter'}
        [m_x] = initParamCALD_x1(vertices, faces);      
    case 'best'
        [neighbors, ~, face_ids] = compute_neigbhors(faces, 2);
        
        [sph_verts] = initParamCALD_x1(vertices, faces);
        [sph_verts_pca] = initParamCALD_x(vertices, faces);
        [cartesian] = start_values_function(vertices, faces, neighbors);
        
        [V, tri] = SpiralSampleSphere(size(vertices, 1));
        hd_gd = hausdorff_distance(V, sph_verts);
        hd_pca = hausdorff_distance(V, sph_verts_pca);
        hd_z = hausdorff_distance(V, cartesian);
        
        [~, min_ind] = min([hd_gd, hd_pca, hd_z]);
        if min_ind == 1
            m_x = sph_verts;
        elseif min_ind == 2
            m_x = sph_verts_pca;
        else
            m_x = cartesian;
        end
    case 'best_reconst'
        [neighbors, ~, face_ids] = compute_neigbhors(faces, 2);
        
        [sph_verts] = initParamCALD_x1(vertices, faces);
        [sph_verts_pca] = initParamCALD_x(vertices, faces);
        [cartesian] = start_values_function(vertices, faces, neighbors);
           
        [is_good, hd_gd] = spherical_parameterization_quality_check(sph_verts, vertices, faces);
        [is_good_1, hd_pca] = spherical_parameterization_quality_check(sph_verts_pca, vertices, faces);
        [is_good_2, hd_z] = spherical_parameterization_quality_check(cartesian, vertices, faces);
     
        [~, min_ind] = min([hd_gd, hd_pca, hd_z]);
        if min_ind == 1
            m_x = sph_verts;
        elseif min_ind == 2
            m_x = sph_verts_pca;
        else
            m_x = cartesian;
        end
    case 'pca'       
        [m_x] = initParamCALD_x(vertices, faces);
    case 'z_axis'
        [neighbors, ~, face_ids] = compute_neigbhors(faces, 2);
        [m_x] = start_values_function(vertices, faces, neighbors);    
end

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

