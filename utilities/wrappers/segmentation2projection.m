function projection = segmentation2projection( filename )

%SEGMENTATION2PROJECTION Creates a projection useful for debugging from the
%segmentation results.

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
%
% Copyright (C) 2014-2016 Murphy Lab
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
% Mar. 4, 2016 I. Cao-Berg If the height of the projection is smaller than
% 100, then the method will not write the filename to the projection

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

if ~isfield( seg, 'cell' )
    warning( 'Segmented cell not present in workspace. Returning empty projection.' );
elseif ~isfield( seg, 'nuc' )
    warning( 'Segmented nucleus not present in workspace. Returning empty projection.' );
else
    c = reshape( seg.cell, size(seg.cell,1), [] );
    d = reshape( seg.nuc, size(seg.nuc,1), [] );
    z = zeros( size(reshape( seg.cell, size(seg.cell,1), [] )) );
    
    position =  [1 1]; % [x y]
    value = [filename];
    projection = cat(3,cat(3,c,d),z);
    
    if size( projection, 1 ) >= 100
        try
            projection = insertText(projection, position, value, ...
                'AnchorPoint', 'LeftTop', ...
                'TextColor', 'white', ...
                'BoxColor', 'black' );
        catch
            disp(['Unable to insert text into projection. ' ...
                'Most likely you are running an older version of Matlab']);
        end
    end
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