function [feat,names] = ml_gmobj2feat(mixes,objintens,param)
%ML_GMOBJ2FEAT Convert Gaussian mixture objects to features.
%   FEAT = ML_GMOBJ2FEAT(MIXES,OBJINTENS) returns a feature vector that is
%   calculated from the gaussian objects MIXES, which is a cell array of
%   Gaussian distributions (See GMM for more details about the data structure).
%   
%   FEAT = ML_GMOBJ2FEAT(MIXES,OBJINTENS,PARAM) specifies how to calculate the
%   features.
%   
%   [FEAT,NAMES] = ML_GMOBJ2FEAT(...) also returns feature names. The features
%   are 'number of objects', 'size parameter', 'intesity parameter'. 
%
%   See also

%   24-Jan-2007 Initial write T. Zhao
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
    error('2 or 3 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

defaultRelation = ml_makecomfun({struct('funname','ml_div2','nvar',2), ...
                   struct('funname','sqrt')});

param = ml_initparam(param,struct('objmodel',struct('stat', ...
    struct('transform',struct('funname','sqrt'),'name','exp'), ...
    'relation',defaultRelation)));

%remove empty structures
rmidx = [];
for i=1:length(mixes)
    if isempty(mixes{i})
        rmidx = [rmidx i];
    end
end

mixes(rmidx) = [];
objintens(rmidx) = [];
feat(1) = length(mixes);
names{1} = 'obj_num';

latents = [];
intensities = [];

for i=1:length(mixes)
    if ~isempty(mixes{i})
        for k=1:mixes{i}.ncentres
            latents = [latents; mixes{i}.covars(1,1,k)];
            intensities = [intensities; ...
               objintens(i)*mixes{i}.priors(k)];
        end
    end        
end

objstat = ml_estpdf(latents,param.objmodel.stat);

%%%%hard-coded%%%%%%%%
feat = [feat objstat.beta];
names{2} = 'objstat_beta';
%%%%%%%%%%%%%%%%%%%%%%%%%

intensdata = ml_evalfun({intensities,latents},param.objmodel.relation);

%%%%hard-coded%%%%%%%%
intenstat = ml_estpdf(intensdata,struct('name','norm'));
feat = [feat intenstat.mu intenstat.sigma];
names{3} = 'intenstat_mu';
names{4} = 'intenstat_sigma';
%%%%%%%%%%%%%%%%%%%%%%%%%%
