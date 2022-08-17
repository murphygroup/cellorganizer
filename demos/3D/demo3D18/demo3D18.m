function answer = demo3D18(options)
% demo3D18
%
% Train 3D generative model of the cell framework (nucleus and cell shape),
% using hole-finding to infer both nucleus and cell shape from the supplied
% protein pattern. The 3D 3T3 dataset was collected in collaboration with
% Dr. Jonathan Jarvik and Dr. Peter Berget.
%
% Input 
% -----
% * a directory of raw or synthetic protein images
% * the resolution of the images (all images should have the same
%   resolution)
%
% Output
% ------
% * a valid SLML model

% Ivan E. Cao-Berg
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

disp( 'demo3D18' );
warning('This demo is deprecated. The demo will be removed in future versions of CellOrganizer');
disp( 'The estimated running time is 2 minutes. Please wait.' );

options.seed = 100;
try
    state = rng( options.seed );
catch err
    rand( 'seed', options.seed ); %#ok<RAND>
end

% generic model options
% ---------------------
% model.name                (optional) Holds the name of the model. Default is empty.
options.model.name = '3D_3T3_framework';

% model.id                  (optional) Holds the id of the model. Default is empty.
options.model.id = num2str(now);

% model.filename            Holds the output filename.
options.model.filename = 'model.xml';

% nuclear shape model options
% ---------------------------
% nucleus.type              Holds the nuclear model type. Default is "medial axis".
options.nucleus.type = 'cylindrical_surface';

% nucleus.id                (optional) Holds the id of the nuclear model. Default is empty.
options.nucleus.id = num2str(now);

% cell shape model options
% ------------------------
% cell.type                 Holds the cell model type. Default is "ratio".
options.cell.type = 'ratio';

% cell.id                   (optional) Holds the id the cell model. Default is empty.
options.cell.id = num2str(now);

% other options
% -------------
% verbose                   (optional) Displays messages to screen. The default is true.
options.verbose = true;

% debug                     (optional) Reports errors and warnings. Default is false.
options.debug = true;

% display                   (optional) Shows displays of intermediate steps
options.display = true;

% train.flag                (optional) Selects what model is going to be trained ('nuclear',
%                           'framework', or 'all'). Default is 'all'
options = ml_initparam( options, struct( ...
    'train', struct( 'flag', 'framework' )));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEEL FREE TO MODIFY THE VARIABLES IN THIS BLOCK

% this demo expects 3D tiff images in this directory and 
% crop files containing a mask for each cell in a masks subdirectory 
imdir = [ pwd filesep '../../../images/3T3/3D/'];
cropdir = [imdir 'masks/'];

files = dir([imdir '*.tif']);
files = sort_nat({files.name});

for i = 1:length(files)
    cell{i} = [imdir files{i}];
end

options.masks = [cropdir filesep '*crop.tif'];

%documentation
options.documentation.description = 'This model has been trained using demo3D18 from CellOrganizer';
options.documentation.dataset = '3D 3T3';
options.documentation.dataset_reference = 'K. Huang and R. F. Murphy (2004) From Quantitative Microscopy to Automated Image Understanding. J. Biomed. Optics 9: 893-912.';
options.documentation.dataset_hyperlink = 'http://murphylab.web.cmu.edu/data/#3D3T3';
options.documentation.author = 'Gregory Johnson';
options.documentation.email = 'gj@andrew.cmu.edu';
options.documentation.copyright = 'Copyright (c) 2007-2017 by the Murphy Lab, Carnegie Mellon University';

%model dimensionality
dimensionality = '3D';

%dataset resolution
options.model.resolution = [0.11, 0.11, 0.5];

%downsampling vector
options.downsampling = [8, 8, 1];

% these control cell segmentation preprocessing (see preprocess)
options.preprocess.stiffness = 0.2; %was originally 0.7
options.preprocess.maxiter = 3000; %same as default
%options.preprocess.quittol = 0.00001; %this is the default
options.preprocess.quittol = 0.0001; %use higher tolerance for speed

% can set this to force single value for nuclear hole finding
% although this is probably not needed (see findDnaMask)
%options.stiffness = ?

% this controls rejecting nuclear hole finding if its volume is too small
options.segminnucfraction = 0.04;

options.nucleus.class = 'nuclear_membrane';
options.nucleus.type = 'cylindrical_surface';

options.cell.class = 'cell_membrane';
options.cell.type = 'ratio';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

answer = img2slml( dimensionality, [], cell, [], options );
end
