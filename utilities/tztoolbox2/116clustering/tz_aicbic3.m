function [aic,bic] = tz_aicbic3(feat,cluster,dist,covtype,ispooled)
%TZ_AICBIC3 AIC BIC calculation with most likely clusters.
%   AIC = TZ_AICBIC3(FEAT,CLUSTER,DIST)
%   
%   AIC = TZ_AICBIC3(FEAT,CLUSTER,DIST,COVTYPE)
%   
%   AIC = TZ_AICBIC3(FEAT,CLUSTER,DIST,COVTYPE,ISPOOLED)
%   
%   [AIC,BIC] = TZ_AICBIC3(...)
%   
%   See also TZ_AICBIC2

%   11-Nov-2006 Initial write T. Zhao
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


if nargin < 3
    error('At least 3 arguments are required');
end

% distance function
if ~exist('dist','var')
    dist = 'eu';
end

%initialize parameters for covariace estimation
initm.method = 'fix';

if ~exist('covtype','var')
    covtype = 0;
end

if ~exist('ispooled','var')
    ispooled = 0;
end

if ispooled~=0 & ispooled~=1
    error('ispooled should be either 0 or 1');
end
initm.ispooled = ispooled;

switch(covtype)
    case 1
        initm.tz_estcov.method = 'mle';
    case 2
        initm.tz_estcov.method = 'ind';
    case 3
        initm.tz_estcov.method = 'idv';     
    case 0
        if strcmp(dist(1:2),'eu')
            initm.tz_estcov.method = 'idv';
        elseif strcmp(dist(1:2),'ma')
            initm.tz_estcov.method = 'mle';
        else
        error(['Wrong type of distance function : ' dist]);
        end
    otherwise
        error(['Unrecognized covariace matrix type: ' num2str(covtype)]);
end

[n, nfeat] = size(feat);

if iscell(cluster)
    clusterid = [];
    for i = 1:length(cluster)
        clusterid(cluster{i})=i;
    end
    ncluster = length(cluster);
    if length(clusterid)~=size(feat,1)
        error(['Number of observations not match in feature matrix and' ...
            ' cluster membership cell array.']);
    end
else
    clusterid = cluster';
    ncluster = max(cluster);
end

initm.labels = clusterid';

%for debugging
% if max(initm.labels==4)
%     keyboard
% end

f = ml_gmmfit(feat,initm);

k = 0;
for i=1:ncluster
    idx = find(clusterid==i);
    if ~isempty(idx)
        k = k+1;
        
        if ~isempty(f.pdfs{k})
            l(k) = ml_loglk(feat(idx,:),f.pdfs{k})+length(idx)*log(f.ps(k)); 
        else
            l(k) = 0;
        end
    end
end 

l = sum(l);

if ~ispooled
    ncov = ncluster;
else
    ncov = 1;
end

p = ncluster-1;

switch initm.tz_estcov.method
    case 'ind'
    % penalty for mean/variance of each cluster
        p = p+ncluster*nfeat+ncov*nfeat;
    case 'mle'
        p = p+ncluster*nfeat+ncov*nfeat*(nfeat+1)/2;
    case 'idv'
        p = p+ncluster*nfeat+ncov;
end

% calculate the IC
aic = 2*p-2*l;
bic = p*log(n)-2*l;
