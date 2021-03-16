function idx = ml_gethridx2(nums,combidx)
%ML_GETHRIDX2 Get herarchical indices.
%   IDX = ML_GETHRIDX2(NUMS,COMBIDX)
%   
%   See also

%   07-Dec-2006 Initial write T. Zhao
%   Copyright (c) 2006 Murphy Lab
%   Carnegie Mellon University
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation; either version 2 of the License,
%   or (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301, USA.
%   
%   For additional information visit http://murphylab.web.cmu.edu or
%   send email to murphy@cmu.edu


if nargin < 2
    error('Exactly 2 arguments are required');
end

offsets = [0 cumsum(nums)];
idxmatrix = repmat(combidx,1,length(offsets));
offsetMatrix = repmat(offsets,length(combidx),1);
diffmatrix = idxmatrix - offsetMatrix;
diffmatrix(diffmatrix<=0) = Inf;
[idx(:,2),idx(:,1)] = min(diffmatrix,[],2);
