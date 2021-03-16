function [v,mu,ps] = ml_estpoolcov(x,label,param)
%ML_ESTPOOLCOV Estimate pooled covariance matrix
%   V = ML_ESTPOOLCOV(X,LABEL)
%   
%   V = ML_ESTPOOLCOV(X,LABEL,PARAM)
%   
%   See also

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


if nargin < 2
    error('2 or 3 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,struct('tz_estcov',struct([])));

removeidx = [];
mu = [];
ps = [];
ncluster = max(label);

for i=1:ncluster
    idx = find(label==i);
    if length(idx)==1
        removeidx = [removeidx;idx];
    end
    if isempty(idx)
        mu(i,:) = zeros(1,size(x,2));
    else
        mu(i,:) = mean(x(idx,:),1);
    end
    
    x(idx,:) = ml_addrow(x(idx,:),-mu(i,:));
    ps(i) = length(idx)/length(label);
end

x(removeidx,:) = [];

v = ml_estcov(x,param.tz_estcov);
 
