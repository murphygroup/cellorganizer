function answer = demo3D33()
% demo3D33
%
% Synthesize multiple 3D images from a lysosome model, at different resolutions.
%
% Input 
% -----
% * a valid CellOrganizer model file
%
% Output
% -------
% * multiple instances of the same cell at different resolutions

% Gregory R. Johnson
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

% July 21, 2012 R.F. Murphy Additional documentation.
% August 10, 2012 I. Cao-Berg Reflected changes from img2projection

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT MODIFY THIS BLOCK
if ~isdeployed()
  current_path = which(mfilename);
  [current_path, filename, extension] = fileparts( current_path );
  cd(current_path);
end

disp( 'demo3D33' );
disp( 'The estimated running time is 30 seconds. Please wait.' );

options.seed = 1245;
options.targetDirectory = pwd;
options.numberOfSynthesizedImages = 1;
options.image.compression = 'lzw';
options.microscope = 'none';
options.sampling.method = 'disc';
options.verbose = true;
options.overlapthresh = 0; %to make the synthesis go very fast
options.prefix = 'imres1';
options.outputres = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    state = rng(options.seed);
catch
    state = RandStream.create('mt19937ar','seed',options.seed);
    RandStream.setGlobalStream(state);
end

if ~exist( [ pwd filesep 'imres1/cell1/cell.tif' ], 'file' )
    slml2img( {'../../../models/3D/lamp2.mat'}, options );
else
    disp('Image exists on disk. Skipping synthesis.');
end

options.prefix = 'imres2';
options.outputres = 0.50;

try
    state = rng(options.seed);
catch
    state = RandStream.create('mt19937ar','seed',options.seed);
    RandStream.setDefaultStream(state);
end
if ~exist( [ pwd filesep 'imres2/cell1/cell.tif'], 'file' )
    answer = slml2img( {'../../../models/3D/lamp2.mat'}, options );
else
    disp('Image exists on disk. Skipping synthesis.');
    answer = 1;
end
