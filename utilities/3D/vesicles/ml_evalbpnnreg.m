function [y,prob] = ml_evalbpnnreg(x,regmodel)

%ML_EVALBPNNREG predicts target values from a BPNN regression model
%   Y=ML_EVALBPNNREG(X,REGMODEL) returns the predicted targets of X.
%   X is a matrix in which rows are samples and columns are variables.
%   Values of Y are in (0,1). REGMODEL is a structure returned from 
%   ML_BPNNREG.
%   If REGMODEL.postp.ctg is 1, it will return labels, which are for 
%   classification.
%   
%   [Y,PROB] = ML_EVALBPNNREG(...) returns a probability matrix.

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

%   14-May-2005 Initial write T. Zhao


if ~strcmp(regmodel.modeltype,'network')
    error('wrong model type: not a network');
end

if isfield(regmodel,'prep')
    if ~isempty(regmodel.prep.featidx)
        x=x(:,regmodel.prep.featidx);
    end
    if regmodel.prep.zscore
        x=ml_zscore(x,regmodel.prep.zmean,regmodel.prep.zsdev);
    end
end

% Summarize network performance
y = mlpfwd(regmodel.trained, x);

%calculate probability
prob = [];
if nargout==2
    
end

if isfield(regmodel,'postp');
    if regmodel.postp.ctg %categorize
        [nmax, y] = max(y,[],2) ;
    end
end
