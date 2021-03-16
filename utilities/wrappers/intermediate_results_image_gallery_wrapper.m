function answer = ...
    intermediate_results_image_gallery_wrapper( input_images_directory, ...
    output_images_directory )

%INTERMEDIATE_RESULTS_IMAGE_GALLERY_WRAPPER Helper method that prepares images
%for an image gallery

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
%
% Copyright (C) 2013-2016 Murphy Lab
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

% Feb. 2, 2016 I. Cao-Berg Updated method to match the refactoring of
% CellOrganizer v2.5
%
% Feb. 2, 2016 I. Cao-Berg Updated method so that it avoids trying to read
% folders as if they were files. The previous iteration worked but displayed
% a lot of warnings. This version is cleaner.
%
% Mar. 4, 2016 I. Cao-Berg Updated method so that it skips synthesis if it
% is unable to generate images from cell parameterization

answer = false;
folders = dir( [ input_images_directory filesep 'param*' ] );

if ~exist( output_images_directory )
    mkdir( output_images_directory )
end

for index=1:1:length(folders)
    disp( ['Processing image file index: ' num2str(index)] );
    try
        folder = folders(index).name;
        folder = [input_images_directory filesep folder ];
        
        if exist( folder, 'dir' )
            [ nuc, cell, improt ] = ...
                intermed2img( folder );
            
            if ~isempty( nuc )
                nuc( find( nuc ~= 0 ) ) = 1;
                cell( find( cell ~= 0 ) ) = 1;
                img = combine( cell, nuc );
                img = im2projection_RGB( img );
                
                value = [ folder ];
                position =  [1 1]; % [x y]
                
                try
                    img = insertText( img , position, value , ...
                        'AnchorPoint', 'LeftTop', ...
                        'TextColor', 'white', ...
                        'BoxColor', 'black' );
                catch
                    if ~exist( 'message', 'var' )
                        message = 'Unable to print text to image. You are probably running an old version of Matlab.';
                        warning( message );
                    end
                end
                
                [ path, filename, ext ] = fileparts( folder );
                filename = ['cell' num2str(sprintf('%05d', str2num(filename(6:end))))];
                output_filename = [ output_images_directory filesep ...
                    filename '.jpg' ];
                imwrite( img, output_filename );
                clear segdna segcell img
            end
        end
    catch err
        disp( 'Unable to synthesize image from cell parameterization.' );
        getReport( err )
    end
end

answer = true;

end%intermediate_results_image_gallery_wrapper

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
