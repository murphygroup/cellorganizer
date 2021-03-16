function [aiclabel,biclabel] = tz_optclusterdir(resultdir,param)
%TZ_OPTCLUSTERDIR Find optimal clusters from a cluster directory.
%   AICLABEL = TZ_OPTCLUSTERDIR(RESULTDIR) returns the [label vector] of the
%   best clustering results based on AICs from the directory RESULTDIR.
%   
%   AICLABEL = TZ_OPTCLUSTERDIR(RESULTDIR,PARAM)
%   
%   [AICLABEL,BICLABEL] = TZ_OPTCLUSTERDIR(...) also returns the [label vector]
%   based on BICs.
%   
%   See also TZ_LOADCLUSTERDIR

%   07-Nov-2006 Initial write T. Zhao
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
    error('1 or 2 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

[aics,bics,nclust] = tz_loadclusterdir(resultdir);

[minaic,k] = min(aics);
k = nclust(k);
tmp = load([resultdir filesep 'cluster' num2str(k)]);
aiclabel = tmp.cluster.label;

[minbic,k] = min(bics);
k = nclust(k);
tmp = load([resultdir filesep 'cluster' num2str(k)]);
biclabel = tmp.cluster.label;
