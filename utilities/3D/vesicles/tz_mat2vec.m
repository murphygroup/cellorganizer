function v=tz_mat2vec(A,bycol)
%TZ_MAT2VEC Convert a matrix to a vector.
%   TZ_MAT2VEC(A) is the same as A(:), which reshape the matrix A to a
%   vector columnwise.
%   
%   TZ_MAT2VEC(A,BYCOL) uses BYCOL to specify columnwise or rowwise
%   convertion. If BYCOL is 0, then the vector will have rowwise order.
%   Otherwise, it is the same as TZ_MAT2VEC(A).

%   ??-???-???? Initial write T. Zhao
%   07-NOV-2004 Modified T. Zhao
%       - add parameters bycol
%   Copyright (c) Murphy Lab, Carnegie Mellon University


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
    error('1 or 2 arguments are required')
end

if ~exist('bycol','var')
    bycol=1;
end

if bycol==0
    A=A';
end

v=A(:);