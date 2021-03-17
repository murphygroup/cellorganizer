function answer = demo2D05(options)
% demo2D05
%
% Train 2D generative pca nuclear and cell shape model using the Murphy Lab 2D HeLa dataset.
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
% Copyright (C) 2018-2019 Murphy Lab
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
  get_murphylab_image_collections( true );
  cd(current_path);
end

disp( 'demo2D05' );
disp( 'The estimated running time is 1 minutes. Please wait...' );
options.verbose = true;
options.debug = false;
options.display = false;
options.save_segmentations = true;

options.model.name = 'demo2D05';
options = ml_initparam( options, struct( ...
    'train', struct( 'flag', 'framework' )));
options.nucleus.id = uuidgen();
options.nucleus.class = 'framework';
options.nucleus.type = 'pca';
options.cell.id = uuidgen();
options.cell.class = 'framework';
options.cell.type = 'pca';

% latent dimension for the model
options.model.pca.latent_dim = 15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dna = {}; cellm = {}; protein = {}; options.masks = {}; options.labels = {};

directory = '../../../images/HeLa/2D/LAM';
for i=1:1:25
    dna{length(dna)+1} = [ directory filesep 'orgdna' filesep 'cell' num2str(i) '.tif' ];
    cellm{length(cellm)+1} = [ directory filesep 'orgcell' filesep 'cell' num2str(i) '.tif' ];
    protein{length(protein)+1} = [ directory filesep 'orgprot' filesep 'cell' num2str(i) '.tif' ];
    options.labels{length(options.labels)+1} = 'LAMP2';
    options.masks{length(options.masks)+1} = [ directory filesep 'crop' filesep 'cell' num2str(i) '.tif' ];
end

directory = '../../../images/HeLa/2D/Nuc';
for i=1:1:25
    dna{length(dna)+1} = [ directory filesep 'orgdna' filesep 'cell' num2str(i) '.tif' ];
    cellm{length(cellm)+1} = [ directory filesep 'orgcell' filesep 'cell' num2str(i) '.tif' ];
    protein{length(protein)+1} = [ directory filesep 'orgprot' filesep 'cell' num2str(i) '.tif' ];
    options.labels{length(options.labels)+1} = 'Nucleoli';
    options.masks{length(options.masks)+1} = [ directory filesep 'crop' filesep 'cell' num2str(i) '.tif' ];
end

options.model.resolution = [ 0.049, 0.049 ];
options.model.filename = 'model.xml';
options.model.id = uuidgen();
options.model.name = '2D HeLa LAMP2+Nuc';
%set nuclei and cell model name
options.nucleus.name = 'LAMP2+Nuc';
options.cell.model = 'LAMP2+Nuc';
%set the dimensionality of the model
dimensionality = '2D';
%documentation
options.documentation.description = 'This model has been trained using demo2D05 from CellOrganizer';
%set model type
options.train.flag = 'framework';

tic; answer = img2slml( dimensionality, dna, cellm, [], options ); toc,
end
