function img = tif2img( filename )
%TIF2IMG Reads a multichannel tif and loads the image to workspace.
%
% Inputs:
% filename = string containing path to image file. files may be of any BioFormats supported type.
%
% Outputs:
% img = resulting image matrix

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
%
% Copyright (C) 2008-2012 Murphy Lab
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

% April 11, 2012 I. Cao-Berg Removed contrast stretch of the image
% June 25, 2012 I. Cao-Berg Updated code so that the Matlab array
%               inherits the properties of the tiff, i.e. the
%               bit-depth

if nargin ~= 1
    error('CellOrganizer: Wrong number of input arguments');
else%loads a multichannel tif to workspace
    
    img = [];
    
    if ~isstr(filename)
        disp('CellOrganizer: filename has to be a string');
        return;
    end
    if ~(exist(filename)==2)
        disp('CellOrganizer: filename does not exist');
        return;
    end
    
    %get information about this image
    info = imfinfo( filename );
    
    %kmliu april 3, 2016 just in case
    if size(info,1)<1
        disp('CellOrganizer: image is empty');
        return;
    end
    
    %kmliu april 3, 2016 remove IMG=[] since this automatically makes its
    %datatype double
    clear img
    %recurse and load every channel to the workspace
    for i=size(info,1):-1:1
        %kmliu april 3, 2016 start from last index to avoid growing array
        %and automatically create matrix with correct datatype
        
        %icaoberg june 25, 2012 by passing the image info
        %the method imread inherits the properties of the original file
        %most importantly the bit depth
        %img(:,:,i) = imread( filename, i ); %#ok<AGROW>
        img(:,:,i) = imread( filename, i, 'Info', info ); %#ok<AGROW>
    end
    
    %icaoberg june 25, 2012
    %don't make the assumption images are uint8
    %img = uint8( img );
    %contrast-stretch the image
    %img = ml_bcimg( img, [], [0 255] );
end
end%tif2img
