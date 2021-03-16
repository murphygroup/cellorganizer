function cofs = ml_gmobjcofs(gaussobjs)
%ML_GMOBJCOFS Find COF of Gaussian objects.
%   
%   See also

%   29-Oct-2006 Initial write T. Zhao
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


if nargin < 1
    error('Exactly 1 argument is required');
end

cofs = {};

for i=1:length(gaussobjs)
    cofs{i} = [];
    for j=1:length(gaussobjs{i})
        if ~isempty(gaussobjs{i}{j})
            for k=1:gaussobjs{i}{j}.ncentres
                cofs{i}(end+1,:) = gaussobjs{i}{j}.centres(k,:);
            end
        end
    end
end

