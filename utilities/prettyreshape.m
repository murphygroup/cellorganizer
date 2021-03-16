function img = prettyreshape( pattern )
% PRETTYRESHAPE Helper method that concatenates multiple reshapes into a
% single image. Useful for displaying a multi channel image into a single
% image. Only works if all images that are identified by the pattern have
% the same size

% Author: Ivan Cao-Berg
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

files = ml_ls( pattern );
if ~isempty( files )
    img = [];
    for i=1:1:length( files )
        img = [ img; im2reshape( ml_readimage( files{i} ) )];
    end
else
    warning(['No files matched the pattern: ' pattern]);
    img = [];
    return
end
end