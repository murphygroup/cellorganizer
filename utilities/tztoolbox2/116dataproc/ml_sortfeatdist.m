function [idx,dists,r] = ml_sortfeatdist(x1,x2)
%ML_SORTFEATDIST Sort features by their distances.
%   IDX = ML_SORTFEATDIST(X1,X2) returns the indices of features with the
%   ascending order of z-scored ascending order between the [feature matrix] 
%   X1 and the [feature matrix] X2. X1 and X2 must have the same number of
%   columns.
%   
%   [IDX,DISTS,R] = ML_SORTFEATDIST(...) also returns the sorted distances
%   and the rank of each feature according to its distances. R(i) is the
%   rank of the ith feature.
%   
%   See also

%   09-Jan-2007 Initial write T. Zhao
%   Copyright (c) 2007 Murphy Lab
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

if isempty(x1) | isempty(x2)
    error('No empty input is allowed');
end

if size(x1,2)~=size(x2,2)
    error('The two input feature matrices must have the same number of columns.')
end

x = zscore([x1;x2]);

x1 = x(1:size(x1,1),:);
x2 = x(size(x1,1)+1:end,:);

mu1 = mean(x1,1);
mu2 = mean(x2,1);

dists = abs(mu1-mu2);
[dists,idx] = sort(dists);

[tmp,r] = sort(idx);
