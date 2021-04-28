function answer = demo3D61()
% demo3D61
%
% Input
% -----
%
% Output
% ------
% * a valid model

% Copyright (C) 2020 Murphy Lab
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT MODIFY THIS BLOCK
if ~isdeployed()
  current_path = which(mfilename);
  [current_path, filename, extension] = fileparts( current_path );
  cd(current_path);
end

disp( 'demo3D61' );
disp( 'The estimated running time is approximately 90 minutes. Please wait.' );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check if dir are present at every run & deletes them to avoid conflicts
% if exist('spharm_result', 'dir')
%     warning('off');
%     rmdir('spharm_result', 's');
%     rmdir('spharm_input', 's');
%     rmdir('param','s');
% end

directory = '../../../images/HeLa/3D/processed';
cell_paths = {}; dna_paths = {}; prot_paths = {}; options.labels = {};
for i = 1:5
    dna_paths{i} = [directory filesep 'LAM_cell' num2str(i) '_ch0_t1.tif'];
    cell_paths{i} = [directory filesep 'LAM_cell' num2str(i) '_ch1_t1.tif'];
    prot_paths{i} = [directory filesep 'LAM_cell' num2str(i) '_ch2_t1.tif'];
    options.labels{length(options.labels)+1} = 'LAMP2';
    options.masks{i} = [directory filesep 'LAM_cell' num2str(i) '_mask_t1.tif'];
end
disp(dna_paths)
pause

%set the dimensionality
dimensionality = '3D';
resolution = [0.049, 0.049, 0.2000];

% parameters for spharm
options_spharm.verbose = true;
options_spharm.debug = ~false;

%set this off as it's causing an error
%Ted 03/18/2020
options_spharm.display = false;

options_spharm = ml_initparam( options_spharm, struct( ...
    'train', struct( 'flag', 'cell' )));
options_spharm.cell.class = 'cell_membrane';
options_spharm.cell.type = 'spharm_rpdm';

% postprocess of parameterization: alignment
options_spharm.model.spharm_rpdm.postprocess = true;
options_spharm.model.resolution = resolution;
options_spharm.downsampling = [1,1,1]; % this is for the spharm modeling of objects
options_spharm.model.filename = 'objects.xml';
options_spharm.model.id = 'objects';
options_spharm.model.name = 'objects';
options_spharm.nucleus.name = 'objects';
options_spharm.cell.model = 'objects';
% degree of the descriptor
options_spharm.model.spharm_rpdm.maxDeg = 31;
% cellular components: either {'cell'}, {'nuc'}, or {'cell', 'nuc'}
options_spharm.model.spharm_rpdm.components = {'cell'};
% latent dimension for the model
options_spharm.model.spharm_rpdm.latent_dim = 15;
% alignment method: 'major_axis' or 'foe'
options_spharm.model.spharm_rpdm.alignment_method = 'major_axis';
% plane of rotation: 'xy' 'yz', 'xz' or 'xyz'
options_spharm.model.spharm_rpdm.rotation_plane = 'xyz';
%documentation
options_spharm.documentation.description = 'This model has been trained for shape-location protein model from CellOrganizer';
options_spharm.model.spharm_rpdm.segminnucfraction = 0.1;
options_spharm.verbose = 1;
options_spharm.spharm_rpdm.NMcost_tol = 1e-7;
options_spharm.spharm_rpdm.NMlargr_tol = 1e-7;
options_spharm.spharm_rpdm.NMfirsttry_maxiter = 300;
options_spharm.spharm_rpdm.NMretry_maxiter = 100;
options_spharm.spharm_rpdm.NMretry_maxiterbig = 300;

% setup_everything:
options.options_spharm = options_spharm;
options.verbose = 0;
%%%
options.train.flag = 'protein';
%options.train.flag = ['framework', 'protein'];
%options.nucleus.class = 'nuclear_membrane';
%options.nucleus.type = 'spharm_rpdm';
%options.cell.class = 'cell_membrane';
%options.cell.type = 'spharm_rpdm';
%%%
options.protein.class = 'vesicle';
options.protein.type = 'spharm_obj'; % new type for vesicle
options.model.id = uuidgen();
options.model.name = 'LAMP2spharm_obj';
options.model.resolution = resolution;
%options.downsampling = [4,4,1];
options.downsampling = [1,1,1]; %this is for the initial input images
options.min_obj_size = 20;
options.max_obj_size = 400;
%options_ppm.min_obj_size=7;
options.local_thresholding_sigma = 5;
options.object_detection_thresPerc = 0.1;
%options.masks = mask_paths;
%options.if_skip_cell_nuclear_model = true;
options.if_skip_cell_nuclear_model = false;

tic; answer = img2slml(dimensionality, dna_paths, cell_paths, prot_paths , options ); toc,

options = [];
options.shape_evolution = 'none';
options.labels = 'unique';
options.subsize = 100; % smaller number means bigger objects
options.viewangle = [0,90]; %down z axis
%options.viewangle = [90,0]; %side view
options.hd_threshold = 10.; % filter out objects with Hausdorff distance greater than this
slml2info({'model.mat'},options);

load('model.mat');
%f = figure('visible','off');
f = figure('visible','on');
hist(log10(model.proteinModel.spharm_obj_model.cellShapeModel.hausdorff_distances),25)
xlabel('Log10 Hausdorff distance between original and SPHARM-RPDM model');
ylabel('Frequency')
saveas( f, 'hausdorff_distance_histogram.png', 'png' );
