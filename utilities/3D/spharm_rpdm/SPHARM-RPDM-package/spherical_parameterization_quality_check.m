function [is_good, hd] = spherical_parameterization_quality_check(sph_verts, vertices, faces, options)
% this function is used to check whether the parameterization is good,
% which could be used to compare different parameterization method. 
% The idea is to perform spharm, and do reconstrucation and see the
% reconstruction error measured by hausdorff distance. 


if nargin < 4
    options = struct();
end

options = process_options_structure(struct('maxDeg', 31, 'hd_thresh', 20));

hd_thresh = options.hd_thresh;

maxDeg = options.maxDeg;
filename = '';
[fvec, deg, Z] = create_SPHARM_des_LSF(vertices, faces, sph_verts, maxDeg, filename);

[vs, fs]=SpiralSampleSphere(size(vertices, 1));
Zs = calculate_SPHARM_basis(vs, deg);
vertices_reconst = real(Zs*fvec);

[hd, lhd, rhd] = hausdorff_distance(vertices, vertices_reconst);

if hd < hd_thresh
    is_good = true;
else
    is_good = false;
end

if false
    figure, patch('vertices', vertices, 'faces', faces, 'FaceVertexCData',jet(size(vertices,1)),'FaceColor','interp');
    view([45, 45]);
    
    [vs, fs]=SpiralSampleSphere(4002);
    Zs = calculate_SPHARM_basis(vs, deg);
    Zvert_pdm = real(Zs*fvec);
    figure, patch('vertices', Zvert_pdm, 'faces', fs, 'FaceVertexCData',jet(size(Zvert_pdm,1)),'FaceColor','interp');
    view([45, 45]);
end

end