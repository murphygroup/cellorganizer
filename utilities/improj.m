function improj( img )
% IMGPROJ Wrapper method for img2projection. Input argument img can
% be a tiff file or a Matlab array.
% 

% Author: Ivan Cao-Berg
% Created: Summer 2012
% 
% Copyright (C) 2012 Murphy Lab
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

if exist( img, 'var' )
    try
    img = img2projection( img );
    imshow( img, [] );
    catch %#ok<*CTCH>
        warning( 'Unable to display projection' );
    end
elseif exist( img, 'file' );
    try
    img = tif2img( img );
    img = img2projection( img );
    imshow( img, [] );
    catch
        warning( 'Unable to open and display tiff file' );
    end
else
    warning('Unrecognized input argument' );
end