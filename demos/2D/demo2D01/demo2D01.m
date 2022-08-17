function answer = demo2D01( options )
% demo2D01
%
% Train a 2D generative model of the nucleus, cell shape, and lysosomes using
% all LAMP2 images in the Murphy Lab 2D HeLa dataset.
%
% Input
% -----
% * a directory of raw or synthetic nucleus images
% * a directory of raw or synthetic cell shape images
% * a directory of raw or synthetic lysosome images
% * the resolution of the images (all images should have the same
%   resolution)
%
% Output
% ------
% * a valid SLML model file

% Ivan E. Cao-Berg (icaoberg@cmu.edu)
%
% Copyright (C) 2013-2020 Murphy Lab
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

disp( 'demo2D01' );
disp( 'The estimated running time is 7-10 minutes. Please wait.' );

options.verbose = true;
options.debug = true;
options.save_segmentations = true;
options.display = false;
options.skip_preprocessing = true;
options.model.name = 'demo2D01';
options = ml_initparam( options, struct( ...
    'train', struct( 'flag', 'all' )));
options.nucleus.class = 'nuclear_membrane';
options.nucleus.type = 'medial_axis';
options.cell.class = 'cell_membrane';
options.cell.type = 'ratio';
options.protein.class = 'vesicle';
options.save_segmentations = true;
options.protein.type = 'gmm';
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

options.temporary_results = [ pwd filesep 'temporary_results' ];
options.model.resolution = [ 0.049, 0.049 ];
options.model.filename = 'lamp2.xml';
options.model.id = 'lamp2';
options.model.name = 'lamp2';

%documentation
options.documentation.description = 'This model has been trained using demo2D01 from CellOrganizer';

answer = img2slml( '2D', dna, cellm, protein, options );
