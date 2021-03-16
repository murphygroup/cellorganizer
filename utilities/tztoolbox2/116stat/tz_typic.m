function [indices,scores] = tz_typic(x,param)
%TZ_TYPIC Select the most typical feature.
%   INDICES = TZ_TYPIC(X) returns the typicalities of the features in the 
%   [feature matrix] X.
%   
%   INDICES = TZ_TYPIC(X,PARAM) speicifies how to calculate the typicalities 
%   by the structure PARAM, which has the following fields:
%       'featsel' - reduce the features or not.
%           'none' : use original data (default)
%           'princomp' : use principal components
%               'numpc' - number of principal components. If numpc is less 
%                   than 1, it means percentage. Default: 0.95.
%               'robust' - robust estimation of the covariance matrix or not.
%                   'true' : robust estimation
%                   'false' : no robust estimation (default)
%
%   [INDICES,SCORES] = TZ_TYPIC(...) also returns the scores.
%   
%   See also

%   28-Feb-2007 Initial write T. Zhao
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

param = ml_initparam(param,struct('featsel','none','numpc',0.95, ...
                                  'robust','false'));

if param.numpc<1
    numpc = -1;
    pcnt = param.numpc;
else
    pcnt = 0;
end

[scores,indices] = ml_TypICfm(x,param.featsel,numpc,pcnt,param.robust);

