function s = ml_sparse( matrix)

% S = ML_SPARSE( MATRIX)
%
% Makes a sparse N-D matrix from a full N-D matrix
%
% Note that matlab's SPARSE/FULL functions only handle 2D matrices.

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

dimensions = size( matrix);
positions = uint32(find( matrix));
values = matrix( positions);
if( length( positions) == 0)
    values = matrix(1);
else
    if( length( find(values == 1)) == length( values))
        values = values(1); % in this case we have a binary image
    end
end
s = struct('dim',dimensions,'pos',positions,'val',values);
