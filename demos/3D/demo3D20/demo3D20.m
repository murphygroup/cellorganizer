function answer = demo3D20(options)
% demo3D20
%
% Train 3D generative diffeomorphic model of the cell framework (nucleus
% and cell shape) using 10 images Murphy Lab 3D HeLa LAMP2 dataset.
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
% -------
% * a valid SLML model file
% * a visualization of the shape space

% Gregory Johnson
%
% Copyright (C) 2013-2017 Murphy Lab
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

disp( 'demo3D20' );
disp( 'The estimated running time is 30 minutes. Please wait.' );

dimensionality = '3D';

% other options
% -------------
options.verbose = true;
options.debug = true;
options = ml_initparam( options, struct( ...
    'train', struct( 'flag', 'framework' )));

% generic model options
% ---------------------
options.model.name = 'demo3D20';
options.model.id = num2str(now);
options.model.filename = 'model.xml' ;

% nuclear shape model options
% ---------------------------
options.nucleus.type = 'diffeomorphic';
options.nucleus.class = 'framework';
options.nucleus.id = num2str(now);
options.nucleus.name = 'lamp2';

% cell shape model options
% ------------------------
options.cell.type = 'diffeomorphic';
options.cell.class = 'framework';
options.cell.model = 'lamp2';
options.cell.id = num2str(now);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEEL FREE TO MODIFY THE VARIABLES IN THIS BLOCK
directory = '../../../images/HeLa/3D/processed';
dna = [directory filesep 'LAM_cell*_ch0_t1.tif'];
cellm = [directory filesep 'LAM_cell*_ch1_t1.tif'];
options.masks = [directory filesep 'LAM_cell*_mask_t1.tif'];
options.model.resolution = [0.049, 0.049, 0.2000];
options.downsampling = [5,5,1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

answer = img2slml( dimensionality, dna, cellm, [], options );
end%demo3D20
