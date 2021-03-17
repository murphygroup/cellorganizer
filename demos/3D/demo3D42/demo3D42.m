function answer = demo3D42()
% demo3D42
%
% This demo illustrates using CellOrganizer to train a protein distribution
% model following the approach described in
%
% K. T. Roybal, T. E. Buck, X. Ruan, B. H. Cho, D. J. Clark, R. Ambler,
% H. M. Tunbridge, J. Zhang, P. Verkade, C. Wülfing, and R. F. Murphy (2016)
% Computational spatiotemporal analysis identifies WAVE2 and Cofilin as 
% joint regulators of costimulation-mediated T cell actin dynamics.  
% Science Signaling 9:rs3. doi: 10.1126/scisignal.aad4149.
%
% The slowest step, which typically takes about 1 min per cell per frame,
% is to align each cell at each time to the standardized template.
% This demo uses 46 cells so it will take about 1 hour on a single core.
%
% Input 
% -----
% * image and annotation files for one or more proteins for one or more
% time points (the default is to use images from the paper of LAT at time 0 
% - downloading the needed images requires about 4 GB of free disk space)
%
% Output
% ------
% * a model for the average concentration in each voxel of a standardized
% cell shape (in demos/LAT_reltime_1.mat)
% * various intermediate results files (in /param and /tmp)

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
  get_murphylab_image_collections( true );
  cd(current_path);
end

% set up intermediate result location
results_location =  './models';

% set up synapse location
synapse_location = {'../../../images/LAT/annotations/*.csv'};      

% options.image_location = image_location;
options.model.tcell.synapse_location = synapse_location;

options.results_location = results_location;
% the running choice for cellorganizer and one sensor of two-point annotation. 
options.model.tcell.named_options_set = 'cellorganizer_one_sensor_2pt';

% set up the mode of synapse to use, as a default, we use one-point, if
% needed you can use two-point by set up the option as true
options.model.tcell.use_two_point_synapses = true;

% set up protein name
options.model.tcell.sensor = 'LAT';

% set up which time points need to analyze
options.model.tcell.timepoints_to_include = [0];

% set up condition 
% options.model.tcell.conditions_to_include = {'Full Stimulus'};
options.model.tcell.model_type_to_include = {'standardized_voxels'};
options.model.resolution = [0.049, 0.049, 0.2000];

options.train.flag = 'protein';

% set up protein model type
options.protein.class = 'standardized_voxels';
options.protein.type = 'standardized_map_half-ellipsoid';

options.train.flag = 'protein';

% set up model name
options.model_prefix = 'LAT_';

options.model.tcell.ometiff = false;
% model.id (optional) Holds the id of the model. Default is empty.

try
    [status, result] = system( 'uuidgen' );
    options.model.id = result(1:end-1);
catch
    options.model.id = num2str(now);
end

% set up model filename
options.model.filename = 'model.xml';

options.debug = true;
options.verbose = true;

dimensionality = '3D';
options.dimensionality = dimensionality;

answer = img2slml(dimensionality, {}, {}, {}, options);
end%demo3D42
