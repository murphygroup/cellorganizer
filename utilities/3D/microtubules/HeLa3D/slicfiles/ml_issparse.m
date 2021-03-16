function status = ml_issparse( matrix)

% STATUS = ML_ISSPARSE( MATRIX)
%
% Tells you whether MATRIX is an MV sparse type or not.

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

if( isfield( matrix, 'dim') & ...
        isfield( matrix, 'pos') & ...
        isfield( matrix, 'val'))
    if( (length(matrix.val) == 1) & (length(matrix.pos)) > 1)
        status = 2; % binary matrix
    else
        status = 1;
    end
else
    status = 0;
end
