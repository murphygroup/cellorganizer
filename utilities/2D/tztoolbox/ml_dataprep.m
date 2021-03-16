function [x,prep] = ml_dataprep(x,param)
%ML_DATAPREP Preprocess data for classification.
%   Y = ML_DATAPREP(X) returns the preprocess data.
%   
%   Y = ML_DATAPREP(X,PARAM) specifies how to preprocess data.
%   
%   [Y,PREP] = ML_DATAPREP(...) also returns the parameters for preprocessing.
%   
%   See also

%   14-Jan-2007 Initial write T. Zhao
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

param = ml_initparam(param,struct('norm',1));

constidx=ml_constidx(x);
if isempty(constidx)
    prep.featidx=[];
else
    allidx=1:size(x,2);
    allidx(constidx)=[];
    prep.featidx=allidx;
    x(:,constidx)=[];
end

prep.zmean = [];
prep.zsdev = [];

%normalize
switch param.norm
    case 1
        [x,zmean,zsdev]=ml_zscore(x);
        prep.zscore=1;
        prep.zmean=zmean;
        prep.zsdev=zsdev;
    case 3
        [x,offset,scale] = ml_libsvmnorm(x);
        prep.zscore = 2;
        prep.zmean = offset;
        prep.zsdev = scale;
    otherwise
        %do nothing
end
