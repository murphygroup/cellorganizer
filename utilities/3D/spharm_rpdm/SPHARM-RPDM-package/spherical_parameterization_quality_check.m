function [is_good, hd] = spherical_parameterization_quality_check(sph_verts, vertices, faces, options)
% this function is used to check whether the parameterization is good,
% which could be used to compare different parameterization method. 
% The idea is to perform spharm, and do reconstrucation and see the
% reconstruction error measured by hausdorff distance. 
% 4/13/2023 R.F. Murphy add figures for debugging
% 5/1/2023 R.F. Murphy don't display debugging figures if deployed
% 5/2/2023 R.F. Murphy fix typo 'ideployed'

if nargin < 4
    options = struct();
end

options = process_options_structure(struct('maxDeg', 31, 'hd_thresh', 20, 'final', false),options);

hd_thresh = options.hd_thresh;

maxDeg = options.maxDeg;
filename = '';
[fvec, deg, Z] = create_SPHARM_des_LSF(vertices, faces, sph_verts, maxDeg, filename);

[vs, fs]=SpiralSampleSphere(size(vertices, 1));
Zs = calculate_SPHARM_basis(vs, deg);
vertices_reconst = real(Zs*fvec);

%figtitle = ['original mesh'];
%figure_filename = sprintf('%soriginal_mesh_%s', figure_dir, imagelegend);
%mesh2figure(vertices,faces,figtitle,figure_filename,dpi);
if options.final && ~isdeployed
    figure(1), patch('vertices', vertices, 'faces', faces, 'FaceVertexCData',jet(size(vs,1)),'FaceColor','interp');
    view([45, 45]);
    figure(2), patch('vertices', vertices_reconst, 'faces', fs, 'FaceVertexCData',jet(size(vs,1)),'FaceColor','interp');
    view([45, 45]);
end

[hd, lhd, rhd] = hausdorff_distance(vertices, vertices_reconst);

if hd < hd_thresh
    is_good = true;
else
    is_good = false;
end

end