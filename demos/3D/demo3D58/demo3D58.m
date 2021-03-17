function answer = demo3D58()
% demo3D58
%
% Synthesize one 3D image with nuclear, cell shape and a vesicular channel using the ratio framework model and write the geometry to a valid VCML file.
%
% Input 
% -----
% * a valid CellOrganizer model file
%
% Output
% ------
% * VCML file
% * three TIFF files (nuclear, cell shape, and lysosome channels)

% Robert F. Murphy
%
% Copyright (C) 2012-2020 Murphy Lab
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
  cd(current_path);
end

start_time = tic;
start_cputime = cputime;
disp( mfilename );
disp( 'The estimated running time is about 2 minutes. Please wait.' );

disp( 'demo3D58' );
options.seed = 12345;
try
    state = rng( options.seed );
catch err
    rand( 'seed', options.seed ); %#ok<RAND>
end

options.targetDirectory = pwd;
options.prefix = 'img';
options.compression = 'lzw';
options.sampling.method = 'disc';
options.debug = false;
options.verbose = true;
options.display = false;
options.output.OMETIFF = true;
options.output.tifimages = true;
options.numberOfGaussianObjects = 20;
options.temporary_results = [ pwd filesep 'temporary_results' ];
options.rendAtStd = 2.0;
options.objstd = options.rendAtStd+0.3;
options.SBML.downsampling = [1/8, 1/8, 1];
options.output.SBMLSpatial = true;
options.SBML.spatial_use_compression = false;
% Branch for this to be fully functional possibly not merged into dev yet but turn on anyway
options.output.VCML.writeVCML = true;
options.VCML.translations = {'cell', 'CP'; 'nucleus', 'NU';  'lamp2_mat', 'EN'; 'CP_EC', 'PM'; 'CP_EN', 'EM'; 'CP_NU', 'NM'};
options.VCML.downsampling = [1/4, 1/4, 1];
options.output.meshes = true;
options.output.indexedimage = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEEL FREE TO MODIFY THE VARIABLES IN THIS BLOCK
list_of_models = {'../../../models/3D/lamp2.mat'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

answer = slml2img( list_of_models, options );


elapsed_time = toc(start_time);
elapsed_cputime = cputime - start_cputime;
fprintf('\n%s took %.3f s (%.3f s CPU time)\n\n', mfilename, elapsed_time, elapsed_cputime);
