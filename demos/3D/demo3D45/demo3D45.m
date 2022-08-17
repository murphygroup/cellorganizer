function answer = demo3D45(options)
% demo3D45
%
% Train 3D generative model of the cell framework (nucleus and cell shape)
% using the Murphy Lab 3D HeLa TfR dataset.
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
% * a valid model

% Xin Lu
%
% Copyright (C) 2012-2017 Murphy Lab
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
%   get_murphylab_image_collections( true );
  cd(current_path);
end

disp( 'demo3D45' );
disp( 'The estimated running time is 20 minutes. Please wait.' );

options.sampling.method = 'disc';
options.debug = true;
options.verbose = true;
options.display = false;

options.downsampling = [5,5,1];
options = ml_initparam( options, struct( ...
    'train', struct( 'flag', 'framework' )));
options.model.filename = 'model.xml';

% generic model options
% ---------------------
options.model.name = '3d_hela_framework_model';
options.model.id = num2str(now);

% nuclear shape model options
% ---------------------------
options.nucleus.type = 'cylindrical_surface';
options.nucleus.class = 'nuclear_membrane';
options.nucleus.name = 'all';
options.nucleus.id = num2str(now);

% cell shape model options
% ------------------------
options.cell.type = 'ratio';
options.cell.class = 'cell_membrane';
options.cell.model = 'framework';
options.cell.id = num2str(now);

options.protein.type = 'gmm';
options.protein.class = 'vesicle';
% options.protein.model = 'framework';
options.protein.id = num2str(now);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEEL FREE TO MODIFY THE VARIABLES IN THIS BLOCK
options.model.resolution = [0.049, 0.049, 0.2000];
imageDirectory = '../../../images/HeLa/3D/processed/';

dimensionality = '3D';
pattern = 'framework';

% dna = [imageDirectory filesep 'TfR*cell*ch0*.tif'];
% cell = [imageDirectory filesep 'TfR*cell*ch1*.tif'];
% options.masks = [imageDirectory filesep 'TfR*mask*.tif'];

directory = '../../../images/ometiff_with_rois';
if ~exist( directory )
    answer = false;
    warning(['Folder ' directory ' does not exist. Exiting demo.']);
    return
end

dna = get_list_of_function_handles_from_wildcards([ directory filesep '/*.ome.tif'], 1);
cell = get_list_of_function_handles_from_wildcards([ directory filesep '/*.ome.tif'], 2);

% documentation
% -------------
options.documentation.author = 'Murphy Lab';
options.documentation.email = 'murphy@cmu.edu';
options.documentation.website = 'murphy@cmu.edu';
options.documentation.description = 'This is the framework model is the result from demo3D45.';
options.documentation.date = date;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

answer = img2slml( dimensionality, dna, cell, [], options );
end%demo3D45