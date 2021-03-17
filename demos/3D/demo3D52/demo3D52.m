function answer = demo3D52(options)
% demo3D52
%
% Train 3D generative SPHARM-RPDM nuclear and cell shape model using the 
% Murphy Lab 3D HeLa dataset.
%
% Input
% -----
% * a directory of raw or synthetic nucleus images
% * a directory of raw or synthetic cell shape images
% * the resolution of the images (all images should have the same
%   resolution)
%
% Output
% ------
% * a valid SLML model file

% Xiongtao Ruan (xruan@andrew.cmu.edu)
%
% Copyright (C) 2018 Murphy Lab
% Computational Biology Department
% School of Computer Science
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu
% 
% 02/24/2019 xruan: specify options for alignment
% 03/05/2021 R.F.Murphy add new options and add displays

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT MODIFY THIS BLOCK
if ~isdeployed()
  current_path = which(mfilename);
  [current_path, filename, extension] = fileparts( current_path );
%  get_murphylab_image_collections( true );
  cd(current_path);
end
disp( 'demo3D52' );
options.verbose = true;
options.debug = false;
options.display = false;
options.model.name = 'demo3D52';
options = ml_initparam( options, struct( ...
    'train', struct( 'flag', 'framework' )));
options.nucleus.class = 'nuclear_membrane';
options.nucleus.type = 'spharm_rpdm';
options.cell.class = 'cell_membrane';
options.cell.type = 'spharm_rpdm';

% criterion for rejecting spherical parameterization
options.hd_thresh = 10;
% postprocess of parameterization: alignment
options.model.spharm_rpdm.postprocess = ~false;
% alignment method: 'major_axis' or 'foe'
options.model.spharm_rpdm.alignment_method = 'major_axis';
% plane of rotation: 'xy' 'yz', 'xz' or 'xyz'
options.model.spharm_rpdm.rotation_plane = 'xyz';

% degree of the descriptor
options.model.spharm_rpdm.maxDeg = 31;
% cellular components: either {'cell'}, {'nuc'}, or {'cell', 'nuc'}
options.model.spharm_rpdm.components = {'cell', 'nuc'};

% latent dimension for the model
options.model.spharm_rpdm.latent_dim = 15;
%increasing these numbers increases compute time but potentially 
% improves model quality; the reported Hausdorff distances between the
% reconstructions and the original shapes can be used to evaluate this
options.NMfirsttry_maxiter = 300;
options.NMretry_maxiter = 100;
options.NMretry_maxiterbig = 300;
%decreasing these numbers decreases compute time but potentially reduces model quality
options.NMcost_tol = 1e-7;
options.NMlargr_tol = 1e-7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the following list of parameters are adapted to the LAMP3 image
% collection, modify these according to your needs
directory = '../../../images/HeLa/3D/processed/';
%dna = [ directory filesep 'LAM_cell[1-9]_ch0_t1.tif' ];
%cellm = [ directory filesep 'LAM_cell[1-9]_ch1_t1.tif' ];
%options.masks = [ directory filesep 'LAM_cell[1-9]_mask_t1.tif' ];
dna = [ directory filesep 'LAM_cell*_ch0_t1.tif' ];
cellm = [ directory filesep 'LAM_cell*_ch1_t1.tif' ];
options.masks = [ directory filesep 'LAM_cell*_mask_t1.tif' ];

options.model.resolution = [0.049, 0.049, 0.2000];
options.downsampling = [5, 5, 1];
options.model.filename = 'lamp2.xml';
options.model.id = 'lamp2';
options.model.name = 'lamp2';
%set nuclei and cell model name
options.nucleus.name = 'LAMP2';
options.cell.model = 'LAMP2';
%set the dimensionality of the model
dimensionality = '3D';
%documentation
options.documentation.description = 'This model has been trained using demo3D52 from CellOrganizer';

options.model.spharm_rpdm.segminnucfraction = 0.1;

tic; answer = img2slml( dimensionality, dna, cellm, [], options ); toc,

options = [];
options.shape_evolution = 'none';
options.labels = 'unique';
options.subsize = 400; % smaller number means bigger objects
%options.viewangle = [0,90]; %down z axis
%options.viewangle = [90,0]; %side view
options.hd_threshold = 10.; % filter out objects with Hausdorff distance greater than this
slml2info({'lamp2.mat'},options);
cd report
web index.html
cd ..

load('lamp2.mat');
f = figure('visible','on');
histogram(log10(model.cellShapeModel.hausdorff_distances(1,:)),25,'FaceColor','red')
hold on
histogram(log10(model.cellShapeModel.hausdorff_distances(2,:)),25,'FaceColor','green')
xlabel('Log10 Hausdorff distance between original and SPHARM-RPDM model');
ylabel('Frequency')
legend('cells','nuclei')
saveas( f, 'hausdorff_distance_histogram.png', 'png' );

end
