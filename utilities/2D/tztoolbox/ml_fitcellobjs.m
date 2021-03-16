function mixes = ml_fitcellobjs(objects,savefile,param)
%ML_FITCELLOBJS Fit objects by gaussian mixtures.
%   
%   See also

%   28-Oct-2006 Initial write T. Zhao
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
    error('2 or 3 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

if exist(savefile,'file')
    tmp = load(savefile);
    mixes = tmp.mixes;
else
    mixes = {};
end

param = ml_initparam(param, ...
    struct('filter',fspecial('gaussian',10,10),'mindist',10, ...
    'isshow',1, 'gmm',struct('covartype','spherical',...
    'options',[0 1 0 0 0 0 0 0 0 0 0 0 0 500])));

if isempty(mixes)
    ncell=1;
    nobj = 0;
else
    ncell = length(mixes);
    nobj = length(mixes{end});
end

for i=ncell:length(objects)   
    for j=nobj+1:length(objects{i})
        [i j]
        obj = objects{i}{j};
        if size(obj,1)==1
            mixes{i}{j} = [];
        else
            mixes{i}{j} = ml_objgaussmix(obj,param);
            if param.isshow==1
                drawnow
            end
            
        end
        save(savefile,'mixes','param');
    end
    nobj = 0;
end
