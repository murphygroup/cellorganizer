function idx = ml_gethridx(nums,combidx)
%ML_GETHRIDX Get hierarchical indices.
%   IDX = ML_GETHRIDX(NUMS,COMBIDX) returns the hierarchial indices of the
%   combined indices COMBIDX according to the cell array NUMS.
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

nums2 = {};
while iscell(nums)   
    [nums,nums2{end+1}] = ml_cellcat(nums);
end
nums2{end+1} = nums;

idx = [];
for i=length(nums2):-1:1
    idx2 = ml_gethridx2(nums2{i},combidx);
    idx = [idx2(:,end) idx];
    combidx = idx2(:,1);
end

idx = [combidx idx];