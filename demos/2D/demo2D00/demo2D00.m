function answer = demo2D00( option )
% demo2D00
%
% Synthesize one 2D image with nuclear, cell shape, and vesicular channels
% from all vesicular object models (nucleoli, lysosomes, endosomes, and
% mitochondria) without convolution. The model was trained from the Murphy
% Lab 2D HeLa dataset.
%
% Input
% -----
% * a list of valid CellOrganizer model files
%
% Output
% ------
% * one TIFF file with six slices (nuclear, cell shape, nucleolar,
%   lysosomal, endosomal, and mitochondrial channels)

% Ivan E. Cao-Berg (icaoberg@cmu.edu)
%
% Copyright (C) 2012-2019 Murphy Lab
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

disp( 'demo2D00' );

options.seed = 12345;
try
    state = rng( options.seed );
catch err
    rand( 'seed', options.seed ); %#ok<RAND>
end

options.targetDirectory = pwd;
options.prefix = 'imgs';
options.compression = 'lzw';
options.debug = false;
options.temporary_results = [ pwd filesep 'temporary_results' ];
options.verbose = false;
options.synthesis = 'all';
options.output.tifimages = false;
options.output.OMETIFF = true;
options.numberOfSynthesizedImages = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

answer = slml2img( {'../../../models/2D/nucleolus.mat', ...
  '../../../models/2D/lysosome.mat', ...
  '../../../models/2D/endosome.mat', ...
  '../../../models/2D/mitochondrion.mat'}, options );
end%demo2D00
