function y = ml_mcfdist(feats,param)
%ML_MCFDIST Pairwise distance between groups of features.
%   Y = ML_MCFDIST(FEATS) returns a vector Y containng the Euclidean distances
%   between each pair of means of feature groups in the MCF FEATS. The meaning
%   of Y is similar to the one returned from PDIST.
%   
%   Y = ML_MCFDIST(FEATS,PARAM) specifies how to calculate the distance by
%   the structure PARAM, which has the following fields:
%       'summary' - summarize each group of features for distance calculation
%           'mean' - mean of features (default).
%           'median' - median of features
%           'mah' - mahalanosis distance
%       'distance' - distance function. This is not available for when
%           'summary' is 'mah'.
%           The choices are from PDIST. See PDIST for more details. The default
%           value is 'euclidean'.
%       'zscored' - zscore the features or not. 1 for zscore and 0 for no 
%           zscore. 
%       'pooled' - this is only used when 'summary' is 'mah'. 1 for pooled 
%           distance and 0 for unpooled distance. The default value is 0.
%
%   See also

%   06-Jan-2007 Initial write T. Zhao
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


if nargin < 1
    error('1 or 2 arguments are required');
end


if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,struct('summary','mean','distance', ...
                                  'euclidean','zscored',1,'pooled',0));

if param.zscored==1
    feats = tz_zscore_mcf(feats);
end

if strcmp(param.summary,'mah')
    nclass = length(feats);
    y = [];
    for i=1:nclass
        for j=i+1:nclass
            y(end+1) = ml_twomaha(feats{i},feats{j},param.pooled,0);
        end
    end
else
    x = [];
    switch param.summary
        case 'mean'
            for i=1:length(feats)
                x(i,:) = mean(feats{i},1);
            end
        case 'median'
            for i=1:length(feats)
                x(i,:) = median(feats{i},1);
            end
        otherwise
            error(['Unrecognized summary option: ' param.summary]);
    end
    y = pdist(x,param.distance);
end


