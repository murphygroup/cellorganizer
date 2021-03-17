function answer = demo3D37(options)
% demo3D37
%
% This demo exists to illustrate how padding size and window size affect the
% performance of diffeomorphic metric.
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

% Copyright (C) 2012-2018 Murphy Lab
% Computational Biology Department
% School of Computer Science
% Carnegie Mellon University

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
answer = false;
if ~isdeployed()
  current_path = which(mfilename);
  [current_path, filename, extension] = fileparts( current_path );
  cd(current_path);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate some superellipses:
 generate_options = struct();
 % Set this:
 generate_options.number_images = 3;
 generate_options.minimum_relative_semidiameter = 1 / 4;
 generate_options.maximum_relative_semidiameter = 2 / 3;
 generate_options.generate_cycle = true;
 % Set this:
 generate_options.image_width = 128;

 [shapes, shape_parameters] = generate_superellipse_shape_set(generate_options);

 options = struct();

% nucleus.id                (optional) Holds the id of the nuclear model. Default is empty.
options.nucleus.type = 'diffeomorphic';
options.nucleus.class = 'framework';
options.cell.type = 'diffeomorphic';
options.cell.class = 'framework';
options.protein.class = 'vesicle';
options.protein.type = 'gmm';

options.verbose = true;
options.debug = true;

% train.flag                (optional) Selects what model is going to be trained ('nuclear',
%                           'framework', or 'all'). Default is 'all'
options = ml_initparam( options, struct( ...
    'train', struct( 'flag', 'framework' )));
options.model.resolution = [1,1,1];

options.model.filename = mfilename;

options.downsampling = [1,1,1];

options.tempparent = [pwd filesep 'temp' filesep];

%documentation
options.documentation.description = 'This model has been trained using demo3D37 from CellOrganizer';

for i = 1:length(shapes)
    objimg{i} = @() shapes{i} > 0;
end

imfunc = @(x) double(shapes{x} > 0);

windowradii = 128./[1,2,4,8,16];

for i = 1:length(windowradii)
    filename = ['diffeomorphic_wr' num2str(windowradii(i))];

    options.model.filename = filename;
    options.model.diffeomorphic.window_radius = windowradii(i);
    options.model.diffeomorphic.image_function = imfunc;

    options.model.diffeomorphic.tempdir = [options.tempparent filesep filename];
    tic
    success = img2slml( '3D', objimg, objimg, [], options );
    timer(i) = toc;
end

answer = true;
end%demo3D37
