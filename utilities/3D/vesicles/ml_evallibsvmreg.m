function [y, prob] = ml_evallibsvmreg(x,regmodel)
%ML_EVALLIBSVMREG classify the input data by a trained LIBSVM model.
%   Y = ML_EVALLIBSVMREG(X,REGMODEL) returns the labels of the [feature matrix]
%   X. These labels are assigned by the trained regression model REGMODEL,
%   which is usually returned from ML_LIBSVMREG.
%   
%   See also ML_LIBSVMREG

%   15-Jan-2007 Initial write T. Zhao
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

if isfield(regmodel,'prep')
    if ~isempty(regmodel.prep.featidx)
        x=x(:,regmodel.prep.featidx);
    end
    if regmodel.prep.zscore
        x=ml_zscore(x,regmodel.prep.zmean,regmodel.prep.zsdev);
    end
end

if nargout==2
    [y, accuracy, prob] = svmpredict(ones(size(x,1),1), x, regmodel.trained, '-b 1');
else
    [y,accuracy] = svmpredict(ones(size(x,1),1), x, regmodel.trained);
end
