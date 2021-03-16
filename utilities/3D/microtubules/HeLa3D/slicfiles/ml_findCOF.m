function COF = ml_findCOF( sparseimg)

% COF = ML_FINDCOF( SPARSEIMG)
%
% Finds the Center of Fluorescence (COF) of the image
% SPARSEIMG. SPARSEIMG must be from ml_sparse.

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

if( ml_issparse( sparseimg))
    % Get graylevel weighted coordinate average
    graylevel = double( sparseimg.val);
    graysum = sum( graylevel);
    [Y,X,Z] = ind2sub( double(sparseimg.dim), double(sparseimg.pos));
    COF(1,1) = sum(Y .* graylevel)/graysum;
    COF(2,1) = sum(X .* graylevel)/graysum;
    COF(3,1) = sum(Z .* graylevel)/graysum;
else
    error('Please input a sparse image!');
end
