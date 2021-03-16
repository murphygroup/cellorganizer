function [y, y_raw] = ml_evalsvmreg(x,regmodel)

%ML_EVALSVMREG classifies the input data by a trained SVM model
%   Y=ML_EVALSVMREG(X,REGMODEL) returns the 1-N coding of X. For example,
%   if there are n classes, Y will be a matrix with n columns. 
%   REGMODEL is a structure returned from ML_SVMREG.
%   If REGMODEL.postp.ctg is 1, it will return labels.
%
%   See also ML_SVMREG
  
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

%HISTORY
%   15-June-2005 Initial write SHANNCC  


if ~strcmp(regmodel.modeltype,'svm')
    error('wrong model type: not a support vector machine');
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
[y, y_raw] = fwd2(regmodel.trained, x);

if isfield(regmodel,'postp');
    if isfield(regmodel.postp,'ctg')
        if regmodel.postp.ctg==1 %categorize
            [nmax, y] = max(y,[],2) ;
        end
    end
end
