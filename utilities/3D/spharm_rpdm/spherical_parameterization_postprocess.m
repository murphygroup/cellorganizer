function [param_postprocess] = spherical_parameterization_postprocess(vertices, faces, sph_verts, fvec, options)
% perform postprocess for spherical parameterization. 
% 
% Author: Xiongtao Ruan
% Date: 09/14/2018
% 
% 02/24/2019 add options for rotation either in xy-plane (default) or xyz 

options = process_options_structure(struct('alignment_method', 'major_axis', ...
                                         'rotation_plane', 'xy', ...    
                                         'maxDeg', 31, ...
                                         'hd_thresh', 50, ...
                                         'use_given_rotation_matrix', false, ...
                                         'use_given_center', false, ...
                                         'check_quality', true), options);

alignment_method = options.alignment_method;
rotation_plane = options.rotation_plane;
maxDeg = options.maxDeg;
check_quality = options.check_quality;
hd_thresh = options.hd_thresh;

[is_good, hd] = spherical_parameterization_quality_check(sph_verts, vertices, faces, options);
if ~is_good
    disp('The parameterization is bad, skip it!');
    param_postprocess = [];
    return;
end

if options.use_given_rotation_matrix
    % assume the bad list for cell and nuclear shapes are the same
    confs.use_given_rotation_matrix = options.use_given_rotation_matrix;
    confs.R = options.rotation_matrix;
end
if options.use_given_center
    confs.use_given_center = options.use_given_center;
    confs.center = options.rotation_center;
end

if strcmp(alignment_method, 'foe') 
    confs.CPoint = 'z';
    % Available values for NPole- 'x';'y';'z'
    confs.NPole = 'x';
    confs.MaxSPHARMDegree = maxDeg;
    confs.rotation_plane = rotation_plane;
    deg = maxDeg;
    [vertices_1, sph_verts_aligned, faces, fvec, fvec_new, R, rotation_center] = align_FOE_2(faces, vertices, sph_verts, fvec, confs);
elseif strcmp(alignment_method, 'major_axis')    
    confs.CPoint = 'z';
    % Available values for NPole- 'x';'y';'z'
    confs.NPole = 'x';
    confs.MaxSPHARMDegree = maxDeg;
    confs.rotation_plane = rotation_plane;    
    deg = maxDeg;
    % [vertices_1, sph_verts_aligned, faces, fvec, fvec_new, R] = align_FOE_2(faces, vertices, sph_verts, fvec, confs);
    [vertices_1, sph_verts_aligned, faces, fvec, fvec_new, R, rotation_center] = align_spharm_by_major_axis(faces, vertices, sph_verts, fvec, confs);            
else
    [sph_verts_aligned] = align_spherical_parameterization(vertices, faces, sph_verts, alignment_method);
    % save(cur_spharm_post_filename, 'sph_verts_aligned');

    surf.vertices=vertices;
    surf.faces=faces;
    sphere.faces = faces;
    sphere.vertices = sph_verts_aligned;
    [surfsmooth, fourier]=SPHARMsmooth2(surf, sphere, maxDeg, 0.001);

    fvec_new = cat(3, fourier.x, fourier.y, fourier.z);
    fvec_new = reshape(fvec_new, [], size(fvec_new, 3));
end

% save(cur_spharm_post_filename, 'sph_verts_aligned', 'fvec_new', 'is_good', 'hd');
% 
% if strcmp(alignment_method, 'foe') || strcmp(alignment_method, 'major_axis')
%     save(cur_spharm_post_filename, 'R', '-append');  
% end

param_postprocess = struct();
param_postprocess.fvec = fvec_new;
param_postprocess.vertices = vertices_1;
param_postprocess.faces = faces;
param_postprocess.sph_verts = sph_verts_aligned;
param_postprocess.hd = hd;
param_postprocess.postprocess = true;

if strcmp(alignment_method, 'foe') || strcmp(alignment_method, 'major_axis')
    param_postprocess.R = R;
    param_postprocess.rotation_center = rotation_center;
end


if false
    [vs fs] = sphereMesh([0 0 0 1]);
    Zs = calculate_SPHARM_basis(vs, deg);

    Zvert = real(Zs*fvec_new);
    figure, patch('vertices', Zvert, 'faces', fs, 'FaceVertexCData',jet(size(vs,1)),'FaceColor','interp');
    view([45, 45]);
end

if false
    [vs, fs]=SpiralSampleSphere(4002);
    Zs = calculate_SPHARM_basis(vs, deg);
    Zvert_pdm = real(Zs*fvec_new);
    figure, patch('vertices', Zvert_pdm(:, [3, 2, 1]), 'faces', fs, 'FaceVertexCData',jet(size(Zvert_pdm,1)),'FaceColor','interp');
    view([45, 45]);
end


end


function [sph_verts_aligned] = align_spherical_parameterization(vertices, faces, sph_verts, alignment_method)
% align mesh by rotation according to the change of north pole. Currently
% implement x-y projection pca
% and just treat the largest z-axis as the north pole

switch alignment_method 
    case 'xy_pca'
        [coeff, score] = pca(vertices(:, 1:2));
        [~, max_ind] = max(score(:, 1));
        new_north_pole = sph_verts(max_ind, :);
    case 'largest_z'
        [~, max_ind] = max(vertices(:, 3));
        new_north_pole = sph_verts(max_ind, :);
end

% Rodrigues' rotation formula
north_pole = [0, 0, 1];
k = cross(north_pole, new_north_pole);

cos_theta = north_pole * new_north_pole';
sin_theta = norm(k);
ku = k ./ sin_theta;
K = [0, -ku(3), ku(2); ku(3), 0, -ku(1); -ku(2), ku(1), 0];
R = eye(3) + sin_theta .* K + (1 - cos_theta) .* K ^ 2;

sph_verts_aligned = sph_verts * R';


end


