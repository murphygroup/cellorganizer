function [n,range] = ml_countnum(x)
%ML_COUNTNUM Count integers in a vector.
%   N = ML_COUNTNUM(X) returns a vector which contains the number of integers
%   occurring in the vector X, which contains integers. n(1) is the number of
%   the minmum in X and n(end) is the number of maxmum in X. 
%   
%   [N,RANGE] = ML_COUNTNUM(X) also returns the range of X.
%    
%   See also

%   28-May-2006 Initial write T. Zhao
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

if nargin < 1
    error('Exactly 1 argument is required')
end

range(1) = min(x);
range(2) = max(x);

n = hist(x,range(2)-range(1)+1);
