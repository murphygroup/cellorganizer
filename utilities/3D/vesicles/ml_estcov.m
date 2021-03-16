function sig = ml_estcov(x,param)
%ML_ESTCOV Covariance matrix estimation.
%   SIG = ML_ESTCOV(X) returns a covariace matrix estimated from the
%   [feature matrix] X.
%   
%   SIG = ML_ESTCOV(X,PARAM) specifies how to calculate the covariance
%   matrix. PARAM has the following fields:
%       'method' - method of estimating covariance matrix
%           'mle' - maximum likelihood estimation
%           'ind' - diagonal covariance matrix. This assumes that all
%               variables are independent.
%           'idv' - This assumes that all variables are independent and
%               have identical variance.
%       'weights' - a colum vector for specifying the weights of the data 
%           points. It must have as many as rows as X. The ith element is
%           the weight for the ith row in X.
%       'calscale' - this is to calibrate the diagonal elements of the
%           covariance matrix is it is not positive definite. The defaut
%           value is 0.01.
%       'mu' - predefined mean of X. If it is empty or not exist, the mean
%           will be calculated from X.
%          
%   See also

%   16-Jun-2006 Initial write T. Zhao
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
    error('1 or 2 arguments are required')
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param, ...
    struct('method','mle','weights',[],'calscale',0.01,'mu',[],'flag',1));
          
if size(x,1)==1
    warning(['There is only one sample for ' ...
        'covariance matrix estimation']);
    sig = 0;
    return;
end
          
if isempty(param.mu)
    if ~isempty(param.weights)
        for i=1:size(x,2)
            param.mu(i) = ml_wmoment(x(:,i),param.weights,1);
        end
    else
        param.mu = mean(x,1);
    end
    
end

switch param.method
    case 'free'
        sig = ml_cov(x,param.flag,param.weights,param.mu);
        [R,p] = chol(sig);
        if p>0 %if the covariance matrix is not positive definite
            warning(['The covariance matrix is not ' ...
                'positive definite. Calibration is done.']);
            sig = sig+diag(diag(sig)*param.calscale);
        end
    case 'mle'
        param.flag = 1;
        param.method = 'free';
        sig = ml_estcov(x,param);

    case 'ind'
        sig = ml_cov(x,param.flag,param.weights,param.mu);
        sig = diag(diag(sig));
%         if isempty(param.weights)
%             sig = diag(var(x,1,1));
%         else
%             for i=1:size(x,2)
%                 sig(i,i) = ml_wmoment(x(:,i),param.weights,2);
%             end
%         end
    case {'idv','iid'}
%         n = size(x,1);
        k = size(x,2);
%         if isempty(param.weights)
%             v = sum(var(x,1,1))/k;
%         else
%             v = ml_wmoment(x(:),repmat(param.weights,k,1),2);
%         end
        x = ml_addrow(x,-param.mu);
        v = ml_cov(x(:),param.flag,repmat(param.weights,k,1),0);
        sig = diag(zeros(1,k)+v);
    otherwise
        error('Unrecognized covariance estimation method.');
end
