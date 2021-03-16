function Distances = ml_eucdist( Center, Points, ScaleFactor)

% DISTANCES = ML_EUCDIST( CENTER, POINTS, SCALEFACTOR)
%
% Finds the euclidian distance of each and every point in POINTS
% (a list of N points, DxN matrix, where columns are [y;x] or 
% [y;x;z] or whatever, depending on D, the number of dimensions),
% with respect to the CENTER point. CENTER is also in the form
% [y;x;z]. SCALEFACTOR is a row vector with D elements, where each
% element specifies a scale factor for scaling the distances in
% each of y, x, z, etc directions. Use this if your pixels are anisotropic.
%
% Created by Meel Velliste in Summer/Fall 2000
% Modified by Meel Velliste in Feb/March 2002 and 6/19/02 to 
% correct for anisotropy

% Copyright (C) 2006  Murphy Lab
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

NPoints = size( Points, 2);
NDims = size( Points, 1);
diff = Points - repmat( Center, [1 NPoints]);

% if( ~exist( 'ScaleFactor', 'var'))
%     if( NDims == 3)
%         ScaleFactor = [1 1 203/48.8];
%     else
%         ScaleFactor = ones(1,NDims);
%     end
% end

diff = diff .* repmat(ScaleFactor',[1 size(diff,2)]);
Distances = sqrt(sum( diff.^2, 1));

