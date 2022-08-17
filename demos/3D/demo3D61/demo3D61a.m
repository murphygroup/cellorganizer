function answer = demo3D61a()
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
% disp( 'The estimated running time is approximately 90 minutes. Please wait.' );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('./ppm'));

cell_paths = {}; dna_paths = {}; prot_paths = {}; options.labels = {};

imageDic = '';
nCells = 100;
directory=strcat(current_path, imageDic);

for i = 0:nCells
    cell_paths{i+1} = [directory filesep sprintf('%03d', i) '_cell.tiff'];
    dna_paths{i+1} = [directory filesep sprintf('%03d', i) '_dna.tiff'];
    prot_paths{i+1} = [directory filesep sprintf('%03d', i) '_protein.tiff'];
    options.labels{length(options.labels)+1} = 'Golgi';
    options.masks{i+1} = [directory filesep sprintf('%03d', i) '_mask.tiff'];
end




%set the dimensionality of the model
dimensionality = '3D';

resolution = [0.29, 0.29, 0.29];



% parameters for spharm
options_spharm.is_demo = false;
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
options_spharm.downsampling = [1,1,1];
options_spharm.model.filename = 'lamp2.xml';
options_spharm.model.id = 'lamp2';
options_spharm.model.name = 'lamp2';
options_spharm.nucleus.name = 'LAMP2';
options_spharm.cell.model = 'LAMP2';
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

% parameters for ppm
options_ppm.datetime_str=strrep(datestr(datetime),' ','_');
options_ppm.aspect=resolution; %shenj's data
options_ppm.min_obj_size=20;
options_ppm.sigma=5;
options_ppm.thresPerc=0.5;
%options_ppm.thresPerc=0.6;
options_ppm.mask_inverted_color_flag=0;
options_ppm.dummy_num='50';
options_ppm.rand_num='70000';
options_ppm.cv_mode='rd_roi';
options_ppm.fold='3';
options_ppm.cv_round='1';
options_ppm.debug=1;

% setup_everything:
options.options_ppm = options_ppm;
options.options_spharm = options_spharm;
options.verbose = 0;
options.train.flag = 'protein';
options.protein.class = 'vesicle';
options.protein.type = 'spharm_obj'; % new type for vesicle
options.model.id = '';
options.model.name = 'spharm_obj';
options.model.resolution = resolution;
options.downsampling = [1,1,1];
options.min_obj_size = 10;
options.max_obj_size = inf;
%options.masks = mask_paths;
options.python_path = '/home/huangqis/miniconda3/bin/python3.9';
options.if_skip_cell_nuclear_model = false;
%options.seg=true;

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
histogram(model.proteinModel.spharm_obj_model.cellShapeModel.hausdorff_distances,25)
set(gca,'YScale','log')
xlabel('Hausdorff distance between original and SPHARM-RPDM model');
ylabel('Frequency')
saveas( f, 'hausdorff_distance_histogram.png', 'png' );
