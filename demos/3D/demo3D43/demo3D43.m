function answer = demo3D43()
% demo3D43
%
% This is the synthesis demo for T cell model. 
% The demo takes in two models: one model contains both cell and nuclear 
% shape models, and the other contains a T cell protein shape model. Same 
% as other synthesis framework, it calls slml2img for the synthesis. The 
% meanings of the options are commented in the script. 
%
% Input 
% -----
% * A protein model with type standardized map halp-elipsoid
% * A framework model the provide the shape of the cell. 
%
% Output
% ------
% * one or more set(s) of synthesized images with cell shape and protein
% pattern. 

% Xiongtao Ruan
%
% Copyright (C) 2016-2017 Murphy Lab
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

if ~isdeployed()
  current_path = which(mfilename);
  [current_path, filename, extension] = fileparts( current_path );
  cd(current_path);
end

% set up intermediate result location
results_location =  './models';

options.results_location = results_location;
% set up the mode of synapse to use, as a default, we use one-point, if
% needed you can use two-point by set up the option as true
options.model.tcell.use_two_point_synapses = true;

% set up protein name
% options.model.tcell.sensor = 'Actin';

% set up which time points need to analyze
options.model.tcell.timepoints_to_include = [0]; 

% set up condition 
% options.model.tcell.conditions_to_include = {'Full Stimulus'};
options.model.tcell.model_type_to_include = {'standardized_voxels'};
options.model.resolution = [0.049, 0.049, 0.2000];

% set up protein model type
options.protein.type = 'morphing';

% set up model name
options.model_prefix = 'LAT_';
% model.id                  (optional) Holds the id of the model. Default is empty.

try
    [status, result] = system( 'uuidgen' );
    options.model.id = result(1:end-1);
catch
    options.model.id = num2str(now);
end

options.model.filename = [ options.model.tcell.model_type_to_include{1}, '.xml' ];

options.debug = false;

dimensionality = '3D';
options.dimensionality = dimensionality;
% cd('../tcell_clean');
current_path = which(mfilename);
[current_path, filename, extension] = fileparts( current_path );
cd(current_path);

% random seeding
options.seed = 12345;
try
    state = rng( options.seed );
catch err
    rand( 'seed', options.seed ); %#ok<RAND>
end

% set up output directory
options.targetDirectory = pwd;

options.prefix = 'img';
% set up compression method for the image
options.compression = 'lzw';
% set up debug
options.debug = false;
% set up verbose
options.verbose = true;
% set up display 
options.display = false;
% set up synthesis method 
options.synthesis = 'all';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up the input models 
% we use a diffeomorphic model for cell and nuclear shape
model_1 = load('../../../models/3D/diffeomorphic/hela_cell_10_15_15.mat');
%model_1 = load('hela_cell_10_15_15.mat');
%T cell model 
try
    model_2 = load('../demo3D42/models/LAT_reltime_0.mat');
catch
    error('You must run demo3D42 first or download a model from the website')
end

% merge two models 
model_2.model.cellShapeModel = model_1.model.cellShapeModel;
model_2.model.nuclearShapeModel = model_1.model.nuclearShapeModel;
model = model_2.model;
% set up model type
model.proteinModel.type = 'standardized_map_half-ellipsoid';
model.proteinModel.class = 'standardized_voxels';
clear model_1
clear model_2

save('model.mat', 'model');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEEL FREE TO MODIFY THE VARIABLES IN THIS BLOCK
list_of_models = {'./model.mat'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output = slml2img( list_of_models, options );

    if (~isempty(output)) 
        answer = true;
    else
        answer = false;
    end
end
