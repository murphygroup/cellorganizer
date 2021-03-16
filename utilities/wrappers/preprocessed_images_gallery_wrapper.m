function answer = ...
    preprocessed_images_gallery_wrapper( input_images_directory, ...
    output_images_directory )

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

% May 17, 2015 I. Cao-Berg Fixed bug in method where it was using output directory
% instead of input directory
%
% Feb. 2, 2016 I. Cao-Berg Updated method to match the refactoring of
% CellOrganizer v2.5
%
% Mar. 4, 2016 I. Cao-Berg Updated method to skip saving if projection is
% already found.

answer = false;
folders = dir( [ input_images_directory filesep 'param*' ] );

if ~exist( output_images_directory )
    mkdir( output_images_directory )
end

for index=1:1:length(folders)
    folder = [ input_images_directory filesep folders(index).name ];
    
    if exist( folder, 'dir' )
        disp( ['Processing image file index: ' num2str(index)] );
        file = [ folder filesep 'seg.mat' ];
        disp( ['Loading image file: ' file] );
        load( file );
        
        method = 1;
        if method == 0
            segdna( find( segdna ~= 0 ) ) = 1;
            segcell( find( segcell ~= 0 ) ) = 1;
            segdna = uint8( segdna );
            segcell = uint8( segcell );
            img = segdna + segcell;
            img = im2projection_RGB( img );
            
            output_filename = [ output_images_directory filesep ...
                folders(index).name '.jpg' ];
            imwrite( img, output_filename );
        else
            projection = segmentation2projection( file );
            output_filename = [ output_images_directory filesep ...
                folders(index).name '.jpg' ];
            if ~isempty( projection )
                if ~exist( output_filename )
                    disp([ 'Saving projection: ' output_filename ]);
                    imwrite( projection, output_filename );
                else
                    disp(['Existing file found: ' output_filename ...
                        ' Skipping computation.'] )
                end
                
            end
        end
    end
end

answer = true;

end%preprocessed_images_gallery_wrapper