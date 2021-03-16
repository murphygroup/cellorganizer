function pts2 = ml_expandpts(pts,order)
%ML_EXPANDPTS Extend points to high order.
%   PTS = ML_EXPANDPTS(PTS,ORDER) returns a matrix, each row of which is
%   extended from a point in PTS by an order of the absolute value of 
%   ORDER. If ORDER is less than 0, zero order will be added.
%   
%   See also

%   07-Apr-2005 Initial write T. Zhao
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

if nargin < 2
    error('Exactly 2 arguments are required')
end

pts2=[];
for i=1:abs(order)
    pts2=[pts2,pts.^i];
end

if order<=0
    pts2=[pts2,ones(size(pts2,1),1)];
end
