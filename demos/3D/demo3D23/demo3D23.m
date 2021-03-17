function answer = demo3D23( options )
% demo3D23
%
% Train 3D generative diffeomorphic nuclear, cell shape, and a
% lysosomal model from all LAMP2 images in the Murphy Lab 3D HeLa dataset.
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

% Gregory Johnson
%
% Copyright (C) 2013-2018 Murphy Lab
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

disp( 'demo3D23' );
warning('This demo is deprecated. The demo will be removed in future versions of CellOrganizer');
disp( 'The estimated running time is 9 hours. Please wait.' );

pattern = 'LAMP2';

directory = '../../../images/HeLa/3D/processed';
dna = [directory filesep 'LAM*cell*ch0*.tif'];
cellm = [directory filesep 'LAM*cell*ch1*.tif'];
protein = [directory filesep 'LAM*cell*ch2*.tif'];
options.masks = [directory filesep 'LAM*cell*mask*.tif'];

dimensionality = '3D';

% options.masks = [directory filesep 'cell*mask*.tif'];

% generic model options
% ---------------------
% model.name                (optional) Holds the name of the model. Default is empty.
options.model.name = 'all';

% model.id                  (optional) Holds the id of the model. Default is empty.
options.model.id = num2str(now);

% model.filename            Holds the output filename.
options.model.filename = [ lower(pattern) '.xml' ];

% nucleus.id                (optional) Holds the id of the nuclear model. Default is empty.
options.nucleus.type = 'diffeomorphic';
options.nucleus.class = 'framework';
options.nucleus.id = num2str(now);

% cell shape model options
% ------------------------
% cell.type                 Holds the cell model type. Default is "ratio".
options.cell.type = 'diffeomorphic';
options.cell.class = 'framework';
options.cell.id = num2str(now);

% protein shape model options
% ---------------------------
% protein.type              (optional) Holds the protein model type. The default is "vesicle".
options.protein.type = 'gmm';
options.protein.class = 'vesicle';
options.protein.id = num2str(now);

% protein.cytonuclearflag   (optional) Determines whether the protein pattern will be generated in
%                           the cytosolic space ('cyto'), nuclear space ('nuc') or everywhere ('all').
%                           Default is cyto.
if strcmpi( options.protein.class, 'nuc' )
    options.protein.cytonuclearflag = 'nuc';
else
    options.protein.cytonuclearflag = 'cyto';
end

% other options
% -------------
% verbose                   (optional) Displays messages to screen. The default is true.
options.verbose = true;

% debug                     (optional) Reports errors and warnings. Default is false.
options.debug = false;

% train.flag                (optional) Selects what model is going to be trained ('nuclear',
%                           'framework', or 'all'). Default is 'all'
options = ml_initparam( options, struct( ...
    'train', struct( 'flag', 'all' )));

%documentation
options.documentation.description = 'This model has been trained using demo3D20 from CellOrganizer';

%model resolution
options.downsampling = [10,10,1];

%this is the resolution of the datasets
options.model.resolution = [0.049, 0.049, 0.2000];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

answer = img2slml( dimensionality, dna, cellm, protein, options );
