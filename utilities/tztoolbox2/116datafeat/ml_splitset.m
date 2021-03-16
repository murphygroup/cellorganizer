function [x1,y1,x2,y2] = ml_splitset(x,y,param)
%ML_SPLITSET
%   X1 = ML_SPLITSET(X,Y)
%   
%   X1 = ML_SPLITSET(X,Y,PARAM)
%   
%   [X1,Y1,X2,Y2] = ML_SPLITSET(...)
%   
%   See also

%   18-Apr-2007 Initial write T. Zhao
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

param = ml_initparam(param,struct('split',.3));

mcf = ml_combfeats2mcf(x,y);
x1 = [];
y1 = [];
x2 = [];
y2 = [];

for i=1:length(mcf)
    n = size(mcf{i},1);
    idx = randperm(n);
    n1 = round(n*param.split);
    n2 = n-n1;
    x1 = [x1;mcf{i}(idx(1:n1),:)];
    y1 = [y1;zeros(n1,1)+i];
    x2 = [x2;mcf{i}(idx(n1+1:end),:)];
    y2 = [y2;zeros(n2,1)+i];
end
