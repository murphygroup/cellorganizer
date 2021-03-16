function data = tz_cleandata(data,rindices)
%TZ_CLEANDATA Clean data.
%   DATA = TZ_CLEANDATA(DATA,RINDICES) removes data points with the indices
%   RINDICES from the matrix or cell array DATA. If DATA is a matrix, the
%   rows with indices in the vector RINDICES will be removed. If DATA is a
%   cell array, RINDICES is also a cell array. RINDICES{I} is a vector and
%   used to clean DATA{I}. If DATA{I} is a cell array, then elements with
%   indices RINDICES{I} will be removed. If DATA{I} is a matrix, rows with
%   indices RINDICES{I} will be removed.
%   
%   See also

%   12-Jul-2006 Initial write T. Zhao
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

if ~iscell(rindices)
    if(size(data,1)==1)
        data(rindices) = [];
    else
        data(rindices,:) = [];
    end
else
    for i=1:length(data)
        if~iscell(data{i})
            data{i}(rindices{i},:) = [];
        else
            data{i}(rindices{i}) = [];
        end
    end
end
