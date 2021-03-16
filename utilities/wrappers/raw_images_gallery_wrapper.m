function answer = ...
    raw_images_gallery_wrapper( input_images_directory, ...
    pattern, output_images_directory, param )
    
%RAW_IMAGES_GALLERY_WRAPPER Helper method that makes useful projections
%from the raw images

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

% Feb. 20, 2014 I. Cao-Berg Inserted try/catch block because method 
% insertText will only work with newer versions of Matlab
%
% Feb. 2, 2016 I. Cao-Berg Updated method to match the refactoring of
% CellOrganizer v2.5

answer = false;

if nargin == 3
    param = [];
end

if ~isfield( param, 'downsample' )
    param = ml_initparam( param, struct( 'downsample', [ 5, 5, 1 ] ) );
end

files =  dir( [ input_images_directory filesep pattern ] );

if ~exist( output_images_directory )
    mkdir( output_images_directory )
end

for index=1:1:length(files)
    output_filename = [ output_images_directory filesep ...
        'cell' num2str( sprintf('%05d', index ) ) '.jpg' ];
    if ~exist( output_filename )
    disp( ['Processing image file index: ' num2str(index)] );
    ch0 = [ input_images_directory filesep files(index).name ];
    ch1 = [ input_images_directory filesep ...
        strrep(files(index).name, 'ch0', 'ch1') ];
    ch2 = [ input_images_directory filesep ...
        strrep(files(index).name, 'ch0', 'ch1') ];
    
    img0 = tif2img( ch0 );
    img0 = ml_downsize( img0, param.downsample, 'linear');
    img0( find( img0 ~= 0 ) ) = 1;
    
    img1 = tif2img( ch1 );
    img1 = ml_downsize( img1, param.downsample, 'linear');
    img1( find( img1 ~= 0 ) ) = 1;

    img2 = tif2img( ch2 );
    img2 = ml_downsize( img2, param.downsample, 'linear');
    img2( find( img2 ~= 0 ) ) = 1;
    
    img = combine3( img1, img2, img0 );
    img = im2projection_RGB( img );
    value = [ 'cell' num2str(index) ...
        '_ch[0,1,2]_t1.tif' ];
    position =  [1 1]; % [x y]
    
    %icaoberg 20/2/2014
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
    
    imwrite( img, output_filename );
    clear img
    else
        disp(['Image ' output_filename ' exists on disk. Skipping processing.'] );
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