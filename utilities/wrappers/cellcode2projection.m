function projection = cellcode2projection( filename )

%CELLCODE2PROJECTION Creates a projection useful for debugging from the
%cell codes temporary file results

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
%
% Feb. 20, 2014 I. Cao-Berg Inserted try/catch block because method insertText will
% only work with newer versions of Matlab

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

if exist( 'cellcodes', 'var' )
    n = [];
    c = [];
    img = [];
    for i=1:1:length(cellcodes)
        n = cat(3,n,cellcodes{i}.nucedge);
        c = cat(3,c,cellcodes{i}.celledge);
        img = [ img, cellcodes{i}.nucedge+cellcodes{i}.celledge ];
    end
    
    n = uint8(n);
    c = uint8(c);
    
    
    img2 = combine( n, c );
    img2 = im2projection_RGB( img );
    
    n = im2projection( n );
    n = imadjust(n,stretchlim(n(:)),[]);
    
    c = im2projection( c );
    c = imadjust(c,stretchlim(c(:)),[]);
    
    projection = [ n, c ];
    projection = cat(3, ...
        cat(3, zeros(size(projection)), projection), ...
        zeros(size(projection)));
    
    position =  [1 1]; % [x y]
    [path, filename, extension] = fileparts( filename );
    filename = [ filename extension ];
    value = [filename];
    
    %icaoberg 20/2/2014
    try
        projection = insertText( projection, position, value, ...
            'AnchorPoint', 'LeftTop', ...
            'TextColor', 'white', ...
            'BoxColor', 'black' );
    catch
        if ~exist( 'message', 'var' )
            message = 'Unable to print text to image. You are probably running an old version of Matlab.';
            warning( message );
        end
    end
else
    warning( 'Variable cellcode not present in workspace. Returning empty projection.')
end

%projection = [ projection, img ];
end%cellcode2projection

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