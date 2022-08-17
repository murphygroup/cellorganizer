function answer = demo3D14()
% demo3D14
%
% Render 2D mean projections along XY, XZ, and YZ axes of images
% synthesized by demo3D00.
%
% Input
% -----
% * a directory of 3D synthetic images
%
% Output
% ------
% * projections of synthetic images as TIFF files

% Ivan E. Cao-Berg
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

disp( 'demo3D14' );
disp( 'The estimated running time is 30 seconds. Please wait.' );

options.verbose = true;
options.debug = true;
options.method = 'mean';
options.compression = 'lzw';
answer=false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEEL FREE TO MODIFY THE VARIABLES IN THIS BLOCK
synthesized_images_directory = '../demo3D00/img/cell1';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist( synthesized_images_directory )
    warning( 'Synthetic images directory does not exist. Exiting demo.' );
    return;
else
    output_folder = [ pwd filesep 'projections' ];
    if ~exist( output_folder )
        mkdir( output_folder );
    end
    
    answer = syn2projection( synthesized_images_directory, ...
        output_folder, options );
end
