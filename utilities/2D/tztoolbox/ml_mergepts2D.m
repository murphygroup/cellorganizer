function pts = ml_mergepts(pts,mindist)
%ML_MERGEPTS Merge close points
%   PTS2 = ML_MERGEPTS(PTS,MINDIST) returns a [point array] base on the
%   input [point array] PTS so that no two points have distance less than
%   MINDIST.
%   
%   See also

%   11-Jul-2006 Initial write T. Zhao
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
    error('Exactly 2 argument is required');
end

dists = squareform(pdist(pts));
dists = dists+diag(zeros(1,size(dists,1))+Inf);
while any(dists(:)<mindist)
    [c,i1,i2] = ml_min(dists);
    pts = [pts;(pts(i1,:)+pts(i2,:))/2];
    pts([i1 i2],:) = [];
    if(size(pts,1)<2)
        break;
    end
    dists = squareform(pdist(pts));
    dists = dists+diag(zeros(1,size(dists,1))+Inf);
end
