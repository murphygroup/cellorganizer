function combidx = ml_getcombidx(nums,idx)
%ML_GETCOMBIDX
%   COMBIDX = ML_GETCOMBIDX(NUMS,IDX)
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

combidx = idx(:,1);

for k=2:size(idx,2)
    if ~isnumeric(nums)
        [nums,nums2] = ml_cellcat(nums);
    else
        nums2 = nums;
    end

    combidx = ml_getcombidx2(nums2,[combidx idx(:,k)]);
end
