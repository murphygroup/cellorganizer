function answer = demo3D21(options)
% demo3D21
%
% Train 3D generative model of the cell framework (nucleus and cell shape),
% using hole-finding to infer both nucleus and cell shape from the supplied
% protein pattern. This is identical to demo3D18 minus scaling the
% images. The 3D 3T3 dataset was collected in collaboration with Dr.
% Jonathan Jarvik and Peter Berget.
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

disp( 'demo3D21' );
warning('This demo is deprecated. The demo will be removed in future versions of CellOrganizer');
disp( 'The estimated running time is 30 minutes. Please wait.' );

%documentation
options.documentation.description = 'This model has been trained using demo3D21 from CellOrganizer';
options.documentation.dataset = '3D 3T3';
options.documentation.dataset_reference = 'K. Huang and R. F. Murphy (2004) From Quantitative Microscopy to Automated Image Understanding. J. Biomed. Optics 9: 893-912.';
options.documentation.dataset_hyperlink = 'http://murphylab.web.cmu.edu/data/3D3T3';
options.documentation.author = 'Gregory Johnson';
options.documentation.email = 'gj@andrew.cmu.edu';
options.documentation.copyright = 'Copyright (c) 2007-2017 by the Murphy Lab, Carnegie Mellon University';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imdir = '../../../images/3T3/3D';
cropdir = [ imdir filesep 'masks' ];
files = dir([imdir filesep '*.tif']);
files = sort_nat({files.name});

for i = 1:length(files)
    cell{i} = @() flipdim(ml_readimage([imdir filesep files{i}]),3);
end

options.masks = [cropdir filesep '*crop.tif'];

% generic model options
% ---------------------
options.model.name = '3D_3T3_framework';
options.model.id = num2str(now);
options.model.filename = 'model.xml';

% nuclear shape model options
% ---------------------------
options.nucleus.type = 'cylindrical_surface';
options.nucleus.class = 'nuclear_membrane';
options.nucleus.id = num2str(now);

% cell shape model options
% ------------------------
options.cell.type = 'ratio';
options.cell.class = 'cell_membrane';
options.cell.id = num2str(now);

% other options
% -------------
options.verbose = true;
options.debug = false;

options = ml_initparam( options, struct( ...
    'train', struct( 'flag', 'framework' )));

%model dimensionality
dimensionality = '3D';

%dataset resolution
options.model.resolution = [0.11, 0.11, 0.5];

%downsampling vector
options.downsampling = [8, 8, 1];

%run demo and make profile
answer = img2slml( dimensionality, [], cell, [], options );
end%demo3D321