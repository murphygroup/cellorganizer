function answer = demo3D32()
% demo3D32
%
% Synthesizes 1 image using a lysosomal model with sampling mode
% set to 'disc', no convolution using the object avoidance methods
% Results will be three TIFF files, one each for cell boundary,
% nuclear boundary, and lysosomes, in folder "synthesizedImages/cell1".
%
% Input
% -----
% * valid SBML file
%
% Output
% ------
% * three TIFF files

% Devin Sullivan
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
  cd(current_path);
end

disp( 'demo3D32' );
disp( 'The estimated running time is 2 minutes 30 seconds. Please wait.' );

answer = false;

options.seed = 3;
try
    state = rng(options.seed);
catch
    state = RandStream.create('mt19937ar','seed',options.seed);
    RandStream.setGlobalStream(state);
end

options.targetDirectory = pwd;
options.prefix = 'img';
options.numberOfSynthesizedImages = 1;
options.compression = 'lzw';
options.microscope = 'none';
options.sampling.method = 'disc';
options.verbose = true;
options.debug = false;
options.output.blenderfile = true;
options.output.blender.downsample = 5;
options.numberOfGaussianObjects = 25;
options.output.tifimages = true;
options.overlapsubsize = 1;
options.overlapthresh = 1;
options.rendAtStd = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

slml2img( {'../../../models/3D/lamp2.mat'}, options );
answer = true;
end%demo3D32
