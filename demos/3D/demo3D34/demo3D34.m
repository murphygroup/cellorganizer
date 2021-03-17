function answer = demo3D34( options )
% demo3D34
%
% Synthesize one 3D image with nuclear, cell shape and a vesicular channel. 
% This demo exports the synthetic image as an OME.TIFF as well as an 
% SBML Spatial instance.
%
% Input 
% -----
% * a valid CellOrganizer model
%
% Output
% ------
% * OME.TIFF
% * SBML instance
% * single channel TIF files

% Ivan E. Cao-Berg
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT MODIFY THIS BLOCK
start_time = tic;
start_cputime = cputime;

if ~isdeployed()
  current_path = which(mfilename);
  [current_path, filename, extension] = fileparts( current_path );
  cd(current_path);
end

disp( 'demo3D34' );
disp( 'The estimated running time is about 31 minutes. Please wait.' );

options.seed = 12345;
try
    state = rng( options.seed );
catch err
    rand( 'seed', options.seed ); %#ok<RAND>
end

options.targetDirectory = pwd;
options.numberOfSynthesizedImages = 1;
options.numberOfGaussianObjects = 25;
options.prefix = 'img';
options.compression = 'lzw';
options.microscope = 'none';
options.sampling.method = 'disc';
options.verbose = true;
options.debug = false;
options.output.tifimages = true;
options.overwrite_synthetic_instances = false;
options.rendAtStd = 1.0;
options.objstd = options.rendAtStd+0.3;
options.overlapsubsize = 1;
options.output.SBMLSpatial = true;
options.output.OMETIFF = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEEL FREE TO MODIFY THE VARIABLES IN THIS BLOCK

% note, the framework will be synthesized using the first protein model
% found to match the given patterns in the SBML file.
% Changing the order/priority of this is not supported at this time.
list_of_models = {'../../../models/3D/tfr.mat'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

answer = slml2img( list_of_models, options );

% validate_SBML_instance is not working
% [sbml_valid, sbml_problems] = validate_SBML_instance('./img/cell1/cell.xml', true)

elapsed_time = toc(start_time);
elapsed_cputime = cputime - start_cputime;
fprintf('\n%s took %.3f s (%.3f s CPU time)\n\n', mfilename, elapsed_time, elapsed_cputime);

end%demo3D34
