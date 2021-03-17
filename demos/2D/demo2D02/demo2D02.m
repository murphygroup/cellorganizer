function answer = demo2D02()
% demo2D02
%
% Synthesize one 2D image with nuclear, cell shape, and lysosomal channels
% from LAMP2 model trained in demo2D01 without convolution.
%
% Input 
% -----
% * a valid CellOrganizer model file
%
% Output
% ------
% * one TIFF file with three slices (nuclear, cell shape, and lysosomal
%   channels)

% Ivan E. Cao-Berg (icaoberg@cmu.edu)
%
% Copyright (C) 2013-2018 Murphy Lab
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

disp( 'demo2D02' );
disp( 'The estimated running time is 1 minute. Please wait.' );

options.seed = 12345;
try
 	state = rng( options.seed );
catch err
 	rand( 'seed', options.seed ); %#ok<RAND>
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%modify the next line to generate more images
options.numberOfSynthesizedImages = 1;

options.targetDirectory = pwd;
options.prefix = 'img';
options.compression = 'lzw';
options.debug = false;
options.verbose = true;
options.display = false;

filename = '../demo2D01/lamp2.mat';
if ~exist( filename )
    warning( [ 'File ' filename ' not found.' ] );
    disp('Using alternative model present in this distribution.');
    filename = '../../../models/2D/endosome.mat';
end

options.output.tifimages = false;
options.output.OMETIFF = true;

slml2img( {filename}, options );
answer = true;
end%demo2D02
