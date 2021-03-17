function answer = demo3D56()
% demo3D56
%
% This demo illustrates using CellOrganizer to train an updated version of
% protein distribution model following the approach described in
%
% K. T. Roybal, T. E. Buck, X. Ruan, B. H. Cho, D. J. Clark, R. Ambler,
% H. M. Tunbridge, J. Zhang, P. Verkade, C. Wülfing, and R. F. Murphy (2016)
% Computational spatiotemporal analysis identifies WAVE2 and Cofilin as
% joint regulators of costimulation-mediated T cell actin dynamics.
% Science Signaling 9:rs3. doi: 10.1126/scisignal.aad4149.
%
% The updates include:
%    1. one point synapse annotation is allowed as valid input;
%    2. a method is implemented for synapse detection with only providing
%       the first time point.
%    3. the method for aligmentment adjustment is implemented.
%
% The slowest step, which typically takes about 1 min per cell per frame,
% is to align each cell at each time to the standardized template.
% This demo uses 46 cells so it will take about 1 hour on a single core.
%
% Input
% -----
% * OMETIFF images with image and annotation files for one or more proteins for the first
% time point (the default is to use images from the paper of LAT at time 0
% - downloading the needed images requires about 4 GB of free disk space)
%
% Output
% ------
% * a model for the average concentration in each voxel of a standardized
% cell shape (in demos/LAT_reltime_1.mat)
% * various intermediate results files (in /param and /tmp)

% Author: Xiongtao Ruan, Xin Lu
%
% Copyright (C) 2016-2019 Murphy Lab
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

% extract tif and annotation files from ometiff
directory = '../../../images/LAT/OME.TIFF/';

% set up intermediate result location
results_location = './models';
options.results_location = results_location;

options.model.tcell.ometiff = true;
% the running choice for cellorganizer and one sensor of two-point annotation.
options.model.tcell.named_options_set = 'cellorganizer_one_sensor_1pt_synapse_infer_aligment_adjust';

% set up the mode of synapse to use, as a default, we use one-point, if
% needed you can use two-point by set up the option as true
options.model.tcell.use_two_point_synapses = ~true;

% set up protein name
options.model.tcell.sensor = 'LAT';

% set up which time points need to analyze
options.model.tcell.timepoints_to_include = [-1, 0, 4];

% set up synapse inference
options.model.tcell.infer_synapses = true;

% set up alignment adjustment
options.model.tcell.adjust_one_point_alignment = true;

% set up condition
% options.model.tcell.conditions_to_include = {'Full Stimulus'};
options.model.tcell.model_type_to_include = {'standardized_voxels'};
options.model.resolution = [0.049, 0.049, 0.2000];

% set up protein model type
options.protein.class = 'standardized_voxels';
options.protein.type = 'standardized_map_half-ellipsoid';

% set up model name
options.model_prefix = 'LAT_';

% model.id (optional)
options.model.id = uuidgen();

% set up model filename
options.model.filename = 'model.xml';

options.debug = true;
options.verbose = true;

dimensionality = '3D';
options.dimensionality = dimensionality;

options.train.flag='protein';

answer = img2slml(dimensionality, {}, {}, {directory}, options);
end
