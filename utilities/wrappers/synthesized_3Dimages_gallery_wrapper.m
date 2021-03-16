function answer = ...
    synthesized_3Dimages_gallery_wrapper( input_images_directory, ...
    output_images_directory )
    
%SYNTHESIZED_3DIMAGES_GALLERY_WRAPPER Helper function that makes colored
%projections from the synthesized images folders used in the demos

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
%
% Copyright (C) 2013 Murphy Lab
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

answer = false;
number_of_files = length( ...
    dir( [ input_images_directory filesep 'cell*' ] ) );

if ~exist( output_images_directory )
    mkdir( output_images_directory )
end

for index=1:1:number_of_files
    disp( ['Processing image file index: ' num2str(index)] );
    current_directory = [ input_images_directory filesep 'cell' num2str(index) ];
    nucleus = [  current_directory ...
        filesep 'nucleus.tif' ];
    cell = [ current_directory ...
        filesep 'cell.tif' ];
    
    img0 = tif2img( nucleus );
    img0( find( img0 ~= 0 ) ) = 1;
    
    img1 = tif2img( cell );
    img1( find( img1 ~= 0 ) ) = 1;
    
    img = combine( img1, img0 );
    img = im2projection_RGB( img );
    
    output_filename = [ output_images_directory filesep ...
        'cell' num2str(index) '_framework.jpg' ];
    imwrite( img, output_filename );
    clear img
    
    counter = 3;
    files = dir( [ current_directory filesep '*.tif' ] );
    for j=1:1:length(files)
        file = files(j).name;
        if strcmpi( file, 'cell.tif' ) || ...
                strcmpi( file, 'nucleus.tif' )
            disp( ['Ignoring ' file ] );
        else
            disp( ['Processing ' file ] );
            img2 = tif2img( [ current_directory filesep ...
                file ] );
            img = combine3( img1, img0, img2 );
            img = im2projection_RGB( img );
            output_filename = [ output_images_directory filesep ...
                'cell' num2str(index) '_ch' num2str(j) '.jpg' ];
            imwrite( img, output_filename );
            clear img
        end
    end
end
answer = true;
end%preprocessed_images_gallery_wrapper

function img = combine( img1, img2 )
img1( find( img1 ~= 0 ) ) = 1;
img2( find( img2 ~= 0 ) ) = 1;

img = img1;
img( find( img2 == 1 ) ) = 2;
end%combine

function img = combine3( img1, img2, img3 )
img1( find( img1 ~= 0 ) ) = 1;
img2( find( img2 ~= 0 ) ) = 1;
img3( find( img3 ~= 0 ) ) = 1;

img = img1;
img( find( img2 == 1 ) ) = 2;
img( find( img3 == 1 ) ) = 3;
end%combine