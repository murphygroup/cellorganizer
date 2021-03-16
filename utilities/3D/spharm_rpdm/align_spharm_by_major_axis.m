function [vertices, sph_verts, faces, fvec, fvec_obj, R, centers_1]=align_spharm_by_major_axis(faces, vertices, sph_verts, fvec, options)
% The function aim for the alignment of both object and parameterization of spharm descriptor, by the major axis of the object in xy plane. 
% The function is copied and modified for align_FOE_3.m 
% 
% Author: Xiongtao Ruan (xruan@cs.cmu.edu)
% Date: 09/17/2018
% 
% Change rotation center as the center in the original space
% 02/24/2019 add options for rotation either in xy-plane (default), xz, yz or xyz 


if nargin < 5
    options = struct();
end

switch deblank(char(options.CPoint))
    case 'x'
        blue = 1;
    case 'y'
        blue = 2;        
    case 'z'
        blue = 3;        
end

switch deblank(char(options.NPole))
    case 'x'
        yellow = 1;
    case 'y'
        yellow = 2;        
    case 'z'
        yellow = 3;        
end

blueyellow = [blue yellow];
degree = options.MaxSPHARMDegree;
rotation_plane = options.rotation_plane;

% rotate in the object space
disp('<< Rotate the object space >>');
% vs = real(Z(:,2:4)*fvec(2:4,:));
% R = object_rotate_R(vs);
% [R, rot_angle] = object_rotate_R_z(vs);
if isfield(options, 'use_given_rotation_matrix') && options.use_given_rotation_matrix && isfield(options, 'R')
    R = options.R;
else
    [R] = object_rotate_pca(vertices, rotation_plane);
end

if isfield(options, 'use_given_center') && options.use_given_center
    center_obj = options.center;
else
    % R = eul2rotm([pi / 4, 0, 0]) * R;
    % R = object_rotate_R([ellipAxes(:,1)';ellipAxes(:,3)']);
    center_obj = real(fvec(1, :)) ./ sqrt(4 * pi);
end

centers_1 = center_obj;
vertices_c = vertices - centers_1;
vertices_c = vertices_c*R'; 
vertices = vertices_c + centers_1;

if false
    % [vs, fs]=SpiralSampleSphere(4002);
    % deg = 31;
    % Zs = calculate_SPHARM_basis(vs, deg);
    % Zvert_pdm = real(Zs*fvec_obj);
    figure, patch('vertices', vertices, 'faces', faces, 'FaceVertexCData',jet(size(vertices,1)),'FaceColor','interp');
    xlabel('x')    
    ylabel('y')
    title('original')
end


% check whether need to flip the meshes
if ~false && ~(isfield(options, 'use_given_rotation_matrix') && options.use_given_rotation_matrix)
    skew = skewness(vertices);

    flipdim = skew(1) < 0;
    if any(flipdim)
        R_skew = diag([-1, -1, 1]);
        % fvec_obj = fvec_obj * R_skew';
        vertices_c = vertices - centers_1;
        vertices_c = vertices_c*R_skew'; 
        vertices = vertices_c + centers_1;  
        R = R_skew * R;
    end
    
    % 02/24/2019 xruan also check the skewness of z
    if strcmp(options.rotation_plane, 'xyz') && skew(3) < 0 
        R_skew = diag([1, -1, -1]);
        % fvec_obj = fvec_obj * R_skew';
        vertices_c = vertices - centers_1;
        vertices_c = vertices_c*R_skew'; 
        vertices = vertices_c + centers_1;  
        R = R_skew * R;        
    end    
end
[fvec_obj, deg, Z_tmp1, name_temp] = create_SPHARM_des_LSF(vertices, [], sph_verts, degree, '', '');        

% rotate the parameter space
disp('<< Rotate the parameter space >>');

vertices_raw = vertices;
sph_verts_raw = sph_verts;
is_upsample = true;
if is_upsample
    % upsample the vertices to improve the performance for finding
    % propriate chart. 
    nvert = size(vertices, 1);
    nvert_upsample = max([min(ceil(nvert * 2 / 1000) * 1000, 2e5), ceil(nvert * 1.5 / 1000) * 1000, min(4e4, ceil(nvert * 5 / 1000) * 1000)]) + 2;
    [vs_upsample, ~]=SpiralSampleSphere(nvert_upsample);
    % deg = d;
    Zs = calculate_SPHARM_basis(vs_upsample, deg);
    Zvert_pdm = real(Zs*fvec_obj);
    sph_verts = vs_upsample;
    vertices = Zvert_pdm;
end

% check if north and south pole is best, almost in xy plane. 
[vertices_pairs] = calc_vertices_pairs(sph_verts);
direction_vecs_raw = vertices(vertices_pairs(:, 1), :) - vertices(vertices_pairs(:, 2), :);
direction_vecs = direction_vecs_raw ./ sqrt(sum(direction_vecs_raw .^ 2, 2));
num_pairs = ceil(size(direction_vecs, 1) * 0.01);
[~, inds] = sort(abs(direction_vecs(:, 3)));
inds_1 = inds(1 : max(num_pairs, sum(abs(direction_vecs(:, 3)) == min(abs(direction_vecs(:, 3)))))); 
[~, max_ind] = max(direction_vecs(inds_1, 1));
% if there are some good ones, then pick the one with longest distance.
x_thrsh = 0.9999;
if sum(direction_vecs(inds_1, 1) > x_thrsh) > 1
    pole_method = 'centroid_closest';
    pole_method = 'longest';
    switch pole_method 
        case 'longest'
            inds_1_max = find(direction_vecs(inds_1, 1) > x_thrsh);
            direction_vec_length = sum(direction_vecs_raw(inds_1(inds_1_max), :) .^ 2, 2);
            [~, max_length_ind] = max(direction_vec_length);
            max_ind = inds_1_max(max_length_ind);
        case 'centroid_closest'
            mean_vertices = mean(vertices);
            inds_1_max = find(direction_vecs(inds_1, 1) > x_thrsh);
            direction_vecs_chosen = direction_vecs(inds_1(inds_1_max), :);
            centroid_vecs = vertices(vertices_pairs(inds_1(inds_1_max), 1), :) - mean_vertices;

            centroid_dists = sum(cross(centroid_vecs, direction_vecs_chosen, 2) .^ 2, 2);

            [~, min_length_ind] = min(centroid_dists);
            max_ind = inds_1_max(min_length_ind);
    end
end
pair_ind = inds_1(max_ind);
pole_inds = vertices_pairs(pair_ind, :);

new_north_pole = sph_verts(pole_inds(1), :);

% Rodrigues' rotation formula
north_pole = [0, 0, 1];
k = cross(north_pole, new_north_pole);
cos_theta = north_pole * new_north_pole';
sin_theta = norm(k);
ku = k ./ sin_theta;
K = [0, -ku(3), ku(2); ku(3), 0, -ku(1); -ku(2), ku(1), 0];
R_param = eye(3) + sin_theta .* K + (1 - cos_theta) .* K ^ 2;
sph_verts = sph_verts * R_param;
sph_verts_raw = sph_verts_raw * R_param;

% The align the equator   
% The idea is find the equator plane, and then find the vertices in the
% same plane as the poles, and find the one with maximum y axis as the
% new equator landmark
pole_vertices = vertices(pole_inds, :);
vert_center_coord = mean(pole_vertices);
% vertices_xy_plane = vertices(:, 3) - mean(pole_vertices(:, 3));
[~, inds] = sort(abs(sph_verts(:, 3)));
inds_1 = inds(1 : num_pairs * 2);
inds_1_pos = inds_1(vertices(inds_1, 2) > vert_center_coord(2));
if isempty(inds_1_pos)
    inds_1_pos = inds_1;
end
vertices_z_diff = abs(vertices(:, 3) - vert_center_coord(3));
[~, inds_2] = sort(vertices_z_diff(inds_1_pos));
num_pairs_1 = max(min(numel(inds_2), 20), sum(vertices_z_diff(inds_1_pos) == min(vertices_z_diff(inds_1_pos))));
inds_3 = inds_2(1 : num_pairs_1);

% make sure the equators are closest to the xy plane of poles, the different of z is minimized. 
% [~, max_ind] = max(vertices(inds_1_pos(inds_3), 2));
[~, min_z_diff_ind] = min(mean(vertices_z_diff(vertices_pairs(inds_1_pos(inds_3), :)), 2));

equator_ind = inds_1_pos(inds_3(min_z_diff_ind));
% check whether the equator index is flipped, which happens in some cases.
if vertices(vertices_pairs(equator_ind, 1), 2) < vertices(vertices_pairs(equator_ind, 2), 2)
    equator_ind = vertices_pairs(equator_ind, 2);
end
new_equator = sph_verts(equator_ind, :);
new_equator_xy = [new_equator(1:2), 0] / norm(new_equator(1:2));
equator = [0, 1, 0];

% Rodrigues' rotation formula
k = cross(equator, new_equator_xy);
cos_theta = equator * new_equator_xy';
sin_theta = norm(k);
ku = k ./ sin_theta;
K = [0, -ku(3), ku(2); ku(3), 0, -ku(1); -ku(2), ku(1), 0];
R_param = eye(3) + sin_theta .* K + (1 - cos_theta) .* K ^ 2;
sph_verts = sph_verts * R_param;
sph_verts_raw = sph_verts_raw * R_param;

% check whether the spherical parameterization is flipped. 
[~, z_ind] = min(pdist2(sph_verts, [1, 0, 0]));
if vertices(z_ind, 3) - vert_center_coord(3) > 0
    sph_verts(:, 1) = -sph_verts(:, 1);
    sph_verts_raw(:, 1) = -sph_verts_raw(:, 1);
end

% [fvec_obj, d, Z_tmp1, name_temp] = create_SPHARM_des_LSF(vertices, [], sph_verts, degree, '', '');        
[fvec_obj, d, Z_tmp1, name_temp] = create_SPHARM_des_LSF(vertices_raw, [], sph_verts_raw, degree, '', '');        

vertices = vertices_raw;
sph_verts = sph_verts_raw;

if false
    [vs, fs]=SpiralSampleSphere(20002);
    % deg = 31;
    Zs = calculate_SPHARM_basis(vs, deg);
    Zvert_pdm = real(Zs*fvec_obj);
    figure, patch('vertices', Zvert_pdm, 'faces', fs, 'FaceVertexCData',jet(size(Zvert_pdm,1)),'FaceColor','interp', 'EdgeColor', [0.1, 0.1, 0.1]);
    xlabel('x')    
    ylabel('y') 
    
    % highlight poles and equator
    points = zeros(6, 3);
    points(1, 3) = 1;  % north pole x+
    points(2, 3) = -1; % south pole x-
    points(3, 2) = -1; % (0, 180) y-
    points(4, 1) = -1; % (0, 270) z+
    points(5, 2) = 1;  % (0, 0)   y+
    points(6, 1) = 1;  % (0, 90)  z-

    dist_mat = pdist2(points, vs);
    [~, inds] = min(dist_mat, [], 2);
    markers = {'*', 'o', '+', 'd', 's', '.'};
    hold on
    for i = 1 : size(points, 1)
        plot3(Zvert_pdm(inds(i), 1), Zvert_pdm(inds(i), 2), Zvert_pdm(inds(i), 3), markers{i}, 'MarkerSize',20, 'MarkerEdgeColor','k', 'LineWidth',3);
        hold on
    end
    hold off;
end


if false
    % sph_verts = sph_verts * diag([1, 1, -1]);
    sph_verts = sph_verts * diag([1, 1, -1]);    
    [fvec_obj_1, d, Z_tmp1, name_temp] = create_SPHARM_des_LSF(vertices, [], sph_verts, degree, '', '');
    [vs, fs]=SpiralSampleSphere(4002);
    deg = 31;
    Zs = calculate_SPHARM_basis(vs, deg);
    Zvert_pdm = real(Zs*fvec_obj_1);
    figure, patch('vertices', Zvert_pdm, 'faces', fs, 'FaceVertexCData',jet(size(Zvert_pdm,1)),'FaceColor','interp');
    xlabel('x')    
    ylabel('y')    
end


end



function [R] = object_rotate_pca(vertices, rotation_plane)
% use PCA method to find the principal axis in xy plane

if nargin < 2
    rotation_plane = 'xy';
end    

switch rotation_plane
    case  'xy'
        inds = [1, 2];
    case 'yz'
        inds = [2, 3];
    case 'xz'
        inds = [1, 3];
    case 'xyz'
        inds = 1 : 3;
end
    
coords = vertices(:, inds);

coeff = pca(coords);
% extremum and saddle points
% expts = U*S;

% 08/30/2018 fix bug for potential reflection in the rotation matrix
if det(coeff) < 0
    if numel(inds) == 2
        coeff = coeff .* [1, -1];
    else
        coeff = coeff .* [-1, -1, -1];
    end
end

% set up rotation matrix
% R = blkdiag(coeff', 1);
R = eye(3);
R(inds, inds) = coeff';

end


function [vertices_pairs] = calc_vertices_pairs(sph_verts)
% make pairs for vertices on the sphere if they form a (approximiate) diameter.
% use knnsearch im matlab instead

n_vert = size(sph_verts, 1);
% vertices_pairs = zeros(n_vert, 2);

% for i = 1 : n_vert
%     vert_i = sph_verts(i, :);
%     [~, pair_i_ind] = min(sum(abs(sph_verts + vert_i), 2));
%     vertices_pairs(i, :) = [i, pair_i_ind];
% end

Idx = knnsearch(sph_verts ,-sph_verts, 'Distance', 'cityblock');

vertices_pairs = [(1 : n_vert)', Idx];

end


