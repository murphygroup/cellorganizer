function projection = segmentation_and_reshape2projection( filename )

%SEGMENTATION_AND_RESHAPE2PROJECTION Creates a projection useful for debugging from the
%segmentation results.

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
%
% Copyright (C) 2014 Murphy Lab
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

projection = [];
if exist( filename )
    try
        load( filename );
    catch
        warning( ['Unable to load file: ' filename '. Returning empty projection' ] );
    end
else
    warning( [ 'File: ' filename ' does not exist. Returning empty projection.' ] );
end

if ~exist( 'segcell', 'var' )
    warning( 'Segmented cell not present in workspace. Returning empty projection.' );
elseif ~exist( 'segdna', 'var' )
    warning( 'Segmented nucleus not present in workspace. Returning empty projection.' );
else
    stacknumber = 15;
    [segdna, segcell] = ml_rescaleImage2Cell( segdna, segcell, stacknumber );
    c = reshape( segcell, size(segcell,1), [] );
    d = reshape( segdna, size(segdna,1), [] );
    z = zeros( size(reshape( segcell, size(segcell,1), [] )) );
    
    position =  [1 1]; % [x y]
    value = [filename];
    projection = cat(3,cat(3,c,d),z);
    %projection = insertText(projection, position, value, 'AnchorPoint', 'LeftTop');
end

end%segmentation2projection

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