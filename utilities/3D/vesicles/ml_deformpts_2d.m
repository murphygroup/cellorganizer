function pts2 = ml_deformpts_2d(pts,A)
%ML_DEFORMPTS_2D Transform 2D points.
%   PTS2 = ML_DEFORMPTS_2D(PTS,A) returns an array of 2d points which is
%   the tranformation of PTS according to transformation matrix A.
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

npt=size(pts,1);

order=-(size(A,1)-1)/2;

qpts=ml_expandpts(pts,order);

pts2=qpts*A;
pts2=pts2(:,1:2);
