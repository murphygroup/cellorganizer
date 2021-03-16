function y = ml_spectral(x,param)
%ML_SPECTRAL data transformation for spectral clustering
%   Y = ML_SPECTRAL(X) returns a [feature matrix] which is the spectral
%   mapping of the input [feature matrix] X.
%   
%   Y = ML_SPECTRAL(X,PARAM) specifies parameters through the fields of
%   PARAM:
%       'scale' - scale for distance. Default value is 1.
%       'k' - number of components to use. Default value is -1, which means
%           all components are used.
%   
%   See also

%   07-Jul-2006 Initial write T. Zhao
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

param = ml_initparam(param,struct('scale',1,'k',-1));

% ds = squareform(pdist(x));

% M = exp(-ds.^2/(2*param.scale^2));

norms = sum(x.^2, 2);
M = repmat(norms, [1 length(norms)]) + repmat(norms', [length(norms)  
1]) - 2 * x * x';
M = exp(-M/(2*param.scale^2));

d = sum(M,1);
sd = d.^-.5;
Q = (sd'*sd).*M;

if param.k<0
    [y,S] = eig(Q);
    y = fliplr(y);
    S = flipud(diag(S));
    y(:,S<0) = [];
    S(S<0) = [];
    param.k = size(y,2);
else
    [y,S] = eigs(Q,param.k);
    S = diag(S);
end

% [y,S] = eig(Q);
% y = fliplr(y);
% S = flipud(diag(S));
% if param.k<0
%     y(:,S<0) = [];
%     S(S<0) = [];
%     param.k = size(y,2);
% else
%     y = y(:,1:param.k);
%     S = y(1:param.k);
% end

y = y.*repmat(sqrt(S)',[size(x,1),1]);

eiglen = sqrt(sum(y.^2,2));
y = y./repmat(eiglen,[1 param.k]);
