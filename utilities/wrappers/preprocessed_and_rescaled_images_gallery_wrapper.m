function answer = preprocessed_and_rescaled_images_gallery_wrapper( ...
    input_images_directory, output_images_directory )

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
files = dir( [ input_images_directory filesep 'cell*.mat' ] );

if ~exist( output_images_directory )
    mkdir( output_images_directory )
end

if ~exist( output_images_directory )
    mkdir( output_images_directory )
end

for index=1:1:length(files)
    file = files(index).name;
    disp( ['Processing image: ' ...
        input_images_directory filesep file] );
    data_filename = [ input_images_directory filesep file ];
   
    [path, filename, extension] = fileparts( data_filename );
    filename = ['cell' num2str(sprintf('%05d', str2num(filename(5:end))))];
    output_filename = [ output_images_directory filesep ...
        filename '.png' ];
    
    if exist( file )
        projection = segmentation2projection( data_filename );
        projection2 = segmentation_and_reshape2projection( data_filename );
        
            padsize = ...
        [0, size(projection,1),abs(size(projection,2)-size(projection2,2))];
    padval = 0;
    
    if size( projection, 2 ) - size( projection2, 2 ) > 0
        %this means projection2 is bigger than projection
        projection3 = padarray( projection2, [0, ...
            abs(size(projection2,2)-size(projection,2))], 'post' );
        projection4 = [ projection; projection3 ];
    elseif size( projection, 2 ) - size( projection2, 2 ) < 0
        %this means projection is bigger than projection2
        projection3 = padarray( projection, [0, ...
            abs(size(projection2,2)-size(projection,2))], 'post' );
        projection4 = [ projection3; projection2 ];
        if ~isempty( projection )
            imwrite( projection4, output_filename );
        end
    else
        %this means they have the same size3
        projection4 = [ projection; projection ];
    end
    
    clear projection projection2 projection3 projection4
end

answer = true;

end%preprocessed_and_rescaled_images_gallery_wrapper