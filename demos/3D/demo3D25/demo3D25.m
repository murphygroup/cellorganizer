function answer = demo3D25()
% demo3D25
%
% Synthesizes 1 image using a lysosomal model with sampling mode
% set to 'disc', no convolution and output.SBML set to true.
% Results will be three TIFF files, one each for cell boundary,
% nuclear boundary, and lysosomes, in folder "synthesizedImages/cell1"
% Additionally, in the folder "synthesizedImages/" will be a
% SBML-Spatial(v0.82a) formatted .xml file containing constructed solid
% geometry(CSG) primitives for lysosomes and parametric objects for the
% cell and nuclear shapes.
%
% These files can then be read into VCell using the built in importer or
% CellBlender using the helper function provided in this distribution.
%
% Input
% -----
% * valid SBML model
%
% Output
% ------
% * three TIFF files
% * XML file with primitives for lysosomes and parametric objects

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
answer = false;
if ~isdeployed()
    current_path = which(mfilename);
    [current_path, filename, extension] = fileparts( current_path );
    cd(current_path);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.seed = 12345;
try
    state = rng( options.seed );
catch err
    rand( 'seed', options.seed ); %#ok<RAND>
end

outputDirectory = pwd;
param = [];
param.targetDirectory = outputDirectory;
param.prefix = 'img';
param.numberOfSynthesizedImages = 1;
param.compression = 'lzw';
param.microscope = 'none';
param.sampling.method = 'disc';
param.verbose = true;
param.debug = false;
param.output.blenderfile = true;
param.output.blender.downsample = 5;
param.output.tifimages = true;
param.output.SBML = true;

answer = slml2img( {'../../../models/3D/lamp2.mat'}, param );