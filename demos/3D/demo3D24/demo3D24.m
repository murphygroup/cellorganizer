function answer = demo3D24( options )
% demo3D24
%
% This demo converts a sample SBML file to an SBML-spatial instance using
% the "matchSBML" function. This function takes an SBML file, matches the
% compartments in the file with available models and synthesizes the
% appropriate instances.
%
% Input
% -----
% * sample SBML file
%
% Output
% ------
% * valid SBML model 

% Devin Sullivan
%
% Copyright (C) 2014-2019 Murphy Lab
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

use_profiling = false;
% use_profiling = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT MODIFY THIS BLOCK
if ~isdeployed()
  current_path = which(mfilename);
  [current_path, filename, extension] = fileparts( current_path );
  cd(current_path);
end

disp( 'demo3D24' );
warning('This demo is deprecated. The demo will be removed in future versions of CellOrganizer');
disp( 'The estimated running time is 1 minute. Please wait.' );

options.seed = 12345;
try
    state = rng( options.seed );
catch err
    rand( 'seed', options.seed ); %#ok<RAND>
end

options.targetDirectory = pwd;
options.numberOfSynthesizedImages = 1;
options.numberOfGaussianObjects = 1;
options.prefix = 'img';
options.compression = 'lzw';
options.microscope = 'none';
options.sampling.method = 'disc';
options.verbose = true;
options.debug = false;
options.output.tifimages = true;
options.rendAtStd = 1;
options.objstd = 1.1;
options.overlapsubsize = 1;

filename = 'Motivating_example_cBNGL2_13_sbml.xml';
url = 'http://www.cellorganizer.org/Downloads/files';

% options.output.SBML = 'Motivating_example_cBNGL2_13_sbml.xml';

options.output.SBMLSpatial = 'Marcon et al. 2016 - Figure 2c Type III.xml';
options.output.SBMLSpatialVCellCompatible = true;
options.output.SBMLSpatialUseCompression = true;
options.output.SBMLSpatialImage = true;
options.output.SBMLDownsampling = [1/16, 1/16, 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEEL FREE TO MODIFY THE VARIABLES IN THIS BLOCK

% note, the framework will be synthesized using the first protein model
% found to match the given patterns in the SBML file.
% Changing the order/priority of this is not supported at this time.
modelpaths = {'../../../models/3D/tfr.mat'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if use_profiling
    profile off; profile('-memory','on'); profile off; profile on;
end

options
options_output = options.output

answer = slml2img(modelpaths,options);

xml_source_filename = [options.prefix, '/cell1/cell.xml'];
xml_destination_filename = [options.prefix, '.xml'];
copyfile(xml_source_filename, xml_destination_filename);

if use_profiling
    profile_results = profile('info');
    profile_filename = ['profile_results', '.mat'];
    profsave(profile_results, [profile_filename, '_html']);
    save(profile_filename, 'profile_results');
end

end%demo3D24