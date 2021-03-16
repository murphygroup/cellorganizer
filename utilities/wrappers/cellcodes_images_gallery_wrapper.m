function answer = ...
    cellcodes_images_gallery_wrapper( input_images_directory, ...
    output_images_directory, param )

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
%
% Copyright (C) 2014 Murphy Lab
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
%
% Feb. 20, 2014 I. Cao-Berg Updated method so it adds trailing zeros to output 
% images filename. This helps for keeping images in order when making automated
% galleries


answer = false;

if nargin == 2
    param = [];
end

if ~isfield( param, 'downsample' )
    param = ml_initparam( param, struct( 'downsample', [ 5, 5, 1 ] ) );
end

files = ...
    dir( [ input_images_directory filesep 'cellcodes_*.mat' ] );

if ~exist( output_images_directory )
    mkdir( output_images_directory )
end

for index=1:1:length(files)
    file = files(index).name;
    disp( ['Processing image: ' ...
        input_images_directory filesep file] );
    data_filename = [ input_images_directory filesep file ];
   
    [path, filename, extension] = fileparts( data_filename );
    %icaoberg 2/20/2014
    filename = [filename(1:10) sprintf('%05d',str2num(filename(11:end)))];
    output_filename = [ output_images_directory filesep ...
        filename '.png' ];
    clear path
    clear filename
    clear extension
    
    projection = cellcode2projection( data_filename );
    if ~isempty( projection )
        imwrite( projection, output_filename );
    end
    clear img
    clear cellcodes
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