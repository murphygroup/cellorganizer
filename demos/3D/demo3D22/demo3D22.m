function answer = demo3D22()
% demo3D22
%
% Synthesizes a protein pattern instance from the synthetic image produced
% in demo3D00.
%
% Input 
% -----
% * a synthetic framework
%
% Output
% ------
% * a synthetic image

% Devin Sullivan
%
% Copyright (C) 2014-2017 Murphy Lab
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

disp( 'demo3D22' );
disp( 'The estimated running time is 30 seconds. Please wait.' );

options.seed = 12345;
try
    state = rng( options.seed );
catch err
    rand( 'seed', options.seed ); %#ok<RAND>
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    folder = '../demo3D00/img';
    copyfile( folder );
catch err
    warning('Unable to copy files. Exiting demo');
    getReport(err)
    return
end

options = [];
options.targetDirectory = pwd;
options.prefix = 'img';
options.numberOfSynthesizedImages = 1;
options.compression = 'lzw';
options.microscope = 'none';
options.sampling.method = 'disc';
options.verbose = true;
options.debug = false;
options.synthesis = 'all';

%make sure objects don't overlap.
options.rendAtStd = 1;
options.objstd = 1.1;
options.overlapsubsize = 1;
options.output.tifimages = 1;
options.numberOfGaussianObjects = 50;

options.model.resolution = [0.049, 0.049, 0.2000];
options.resolution.cell = [0.049, 0.049, 0.2000];
img = tif2img( [ pwd filesep 'cell1' filesep 'nucleus.tif' ] );
options.instance.nucleus = img;
options.instance.resolution= [0.049, 0.049, 0.2000];
clear img

img = tif2img( [ pwd filesep 'cell1' filesep 'cell.tif' ] );
options.instance.cell = img;
clear img

answer = slml2img( {'../../../models/3D/tfr.mat'}, options );
end%demo3D22
