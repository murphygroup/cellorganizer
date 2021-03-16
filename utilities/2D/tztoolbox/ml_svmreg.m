function regmodel = ml_svmreg(x,y,t)
%ML_BPNNREG Multivariate regression by SVM
%   REGMODEL=ML_BPNNREG(X,Y,T) returns a structure restoring 
%   the regression model of f(X)=E(Y|X). T is the structure of
%   parameters. It has the following fields:
%       norm, stop, randtrainsel: see ML_REGRESS
%       model_types, tutor,C_values, rbf_levels: see SVM toolbox
%   Y could be one column or n column. For example,
%           Y=[1 1 2 2 3 3]';
%           is the same as
%           Y=[1 -1 -1;1 -1 -1;-1 1 -1;-1 1 -1;-1 -1 1;-1 -1 1];
%
%   SEE ALSO ML_EVALSVMREG

% Copyright (C) 2006  Murphy Lab
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

%   15-June-2005 written by SHANNCC

if ~exist('t','var')
    t.norm=1;
end

if ~isfield(t,'norm')
    t.norm=1;
end

if ~isfield(t,'stop')
    t.stop=0;
end

prep.zmean=[];
prep.zsdev=[];
prep.zscore=0;
constidx=ml_constidx(x);
if isempty(constidx)
    prep.featidx=[];
else
    allidx=1:size(x,2);
    allidx(constidx)=[];
    prep.featidx=allidx;
    x(:,constidx)=[];
end

%normalize upon all training data
if t.norm==1
    [x,zmean,zsdev]=ml_zscore(x);
    prep.zscore=1;
    prep.zmean=zmean;
    prep.zsdev=zsdev;
end

if t.stop==1

else
    net = t.model_types;
    net = train(t.model_types, t.tutor, x, y, t.C_values, rbf(t.rbf_levels));
end

regmodel=struct('modelname','svm','modeltype',...
    'svm','trained',net,'t',t,'prep',prep);

if size(y,2)==1
    regmodel.postp.ctg=1;
end
