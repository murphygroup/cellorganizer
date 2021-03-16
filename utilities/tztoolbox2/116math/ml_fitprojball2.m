function f = ml_fitprojball2(x,y,param)
%ML_FITPROJBALL2
%   PARA = ML_FITPROJBALL2(X,Y)
%   
%   See also

%   29-Mar-2007 Initial write T. Zhao
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
    error('Exactly 2 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,struct('minv',[1 1]));

para = lsqcurvefit(@fitfun,[5 1.5], ...
    x,y,param.minv,[Inf,Inf]);
para = lsqcurvefit(@fitfun,para, ...
    x,y,param.minv,[Inf,Inf]);
f = struct('funname','ml_projball2','param',para(2));
f.transform.output.param.scale = para(1);
f.transform.output.funname = 'ml_linfun';

%%%%%%%%%%%%%%%%%%%
function y = fitfun(t,x)

y = t(1)*ml_evalfun(x,struct('funname','ml_projball2','param',t(2)));
