function [ncm,pcm] = ml_heteroconfmat(y1,y2)
%ML_HETEROCONFMAT Hetero-confusion matrix.
%   NCM = ML_HETEROCONFMAT(Y1,Y2) returns a table of integers. The element at
%   the ith row and jth column is the number of occurrence with label i in
%   the [label vector] Y1 and label j in the [label vector] Y2. So the size of
%   table, or confusion matrix, is MxN, where M is the maximum label in Y1 and 
%   N is the maximum label in Y2.
%   
%   [NCM,PCM] = ML_HETEROCONFMAT(...) also returns the confusion matrix with
%   percentages.
%   
%   See also

%   05-Dec-2006 Initial write T. Zhao
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

nclass1 = max(y1);
nclass2 = max(y2);

ncm = zeros(nclass1,nclass2);

for i=1:nclass1
    idx = find(y1==i);
    if ~isempty(idx)
        [n,range] = ml_countnum(y2(idx));
        ncm(i,range(1):range(2)) = n;
    end
end

if nargout>=2
    pcm = ml_normrow(ncm)*100;
end


