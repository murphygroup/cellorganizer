function answer = demo3D63( options )
% demo3D63
%
% Synthesize one 3D image with nuclear and cell shape channels using the
% SPHARM framework model and combine with reaction network in VCML format
% into a single VCML file.
%
% Input
% -----
% * a valid CellOrganizer model
% * one or more valid VCML reaction network models with Virtual Cell's
%   default units
%
% Output
% ------
% * VCML file
% * single channel TIF files

% Author: Ivan E. Cao-Berg, Taraz Buck
%
% Copyright (C) 2016-2020 Murphy Lab
% Lane Center for Computational Biology
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
% FEEL FREE TO MODIFY THE VARIABLES IN THIS BLOCK

% The framework (cell and nucleus shapes) will be synthesized using framework_model
framework_model = '../../../models/3D/spharm/lamp2.mat';
reaction_network_file = '../../../data/SarmaGhosh2012ForCOdiff1e-2.vcml';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT MODIFY THIS BLOCK
if ~isdeployed()
    current_path = which(mfilename);
    [current_path, filename, extension] = fileparts( current_path );
    cd(current_path);
end

start_time = tic;
start_cputime = cputime;
disp( mfilename );
disp( 'The estimated running time is about 3 minutes. Please wait.' );

options.seed = 741332;
options.targetDirectory = pwd;
options.prefix = 'img';
options.synthesis = 'framework';
options.model.spharm_rpdm.synthesis_method = 'random_sampling';
options.model.spharm_rpdm.imageSize = [512, 512, 36];
options.model.spharm_rpdm.synthesis_resolution = [0.10 0.10 0.1];
options.numberOfSynthesizedImages = 1;

options.compression = 'lzw';
options.microscope = 'none';
options.sampling.method = 'disc';
options.verbose = true;
options.debug = false;
options.output.tifimages = true;
options.output.OMETIFF = true;
options.output.indexedimage = true;

% VCML Options
options.output.VCML.writeVCML = true;
options.VCML.downsampling = 1/2;
options.VCML.translations = {'cell', 'CP'; 'nuc', 'NU'; ...
    'nucleus', 'NU'; 'lamp2_mat_tfr_mat', 'EN'; 'CP_EC', 'PM'; ...
    'CP_EN', 'EM'; 'CP_NU', 'NM'};

options.VCML.input_filename = reaction_network_file;
options.NET.useImageAdjacency = false;

options.VCML.default_time_step = 1e0;
options.VCML.end_time = 4000;
options.VCML.output_time_step = 100;
options.VCML.max_time_step = 4;
options.VCML.absolute_tolerance = 1e-8;
options.VCML.relative_tolerance = 1e-8;

try
    state = rng( options.seed );
catch err
    rand( 'seed', options.seed ); %#ok<RAND>
end

% Synthesize Image and VCML file
answer = slml2img( { framework_model }, options );
end
