function B = ml_multrow(A,x)
%ML_MULTROW Multiply the a row vector to each row of a matrix.
%   ML_MULTROW(A,X) returns a matrix in which each row is the array product
%   between the corresponding row in matrix A and the row vector X. A and X
%   must have the same number of columns if X is not a scalar.
%   
%   See also ML_ADDROW

%   13-Mar-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

% Copyright (C) 2007  Murphy Lab
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

if nargin < 2
    error('Exactly 2 arguments are required')
end

if length(x(:))==1
    B = A*x;
else
    if(size(A,2)~=size(x,2))
        error('The matrix and the vector must have the same number of columns');
    end
    B=A.*(ones(size(A,1),1)*x);
end