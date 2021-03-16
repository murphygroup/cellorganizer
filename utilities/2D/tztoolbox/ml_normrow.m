function p = ml_normrow(x)
%ML_NORMROW Normalize data so that sum of each row is 1.
%   P = ML_NORMROW(X) returns a matrix that has the same size as X. The
%   sum of each row in P is 1. Basically, P(i,j) = X(i,j)/sum(X(i,:)).
%   If the sum of row i is 0, then the ith of P is also a zero vector.
%
%   See also

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

%   05-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

rowSum = sum(x,2);
nonzeroIndices = find(rowSum~=0);
p = zeros(size(x));
rowSum = rowSum(nonzeroIndices);

p(nonzeroIndices,:) = x(nonzeroIndices,:)./ ...
    repmat(rowSum,1,size(x,2));