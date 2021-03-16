function sig = tz_estcov(x,param)
%TZ_ESTCOV Obsolete. See ML_ESTCOV.
%   SIG = TZ_ESTCOV(X) returns a covariace matrix estimated from the
%   [feature matrix] X.
%   
%   SIG = TZ_ESTCOV(X,PARAM) specifies how to calculate the covariance
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

error(tz_genmsg('of','tz_estcov','ml_estcov'));

if nargin < 1
    error('1 or 2 arguments are required')
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param, ...
    struct('method','mle','weights',[],'calscale',0.01));
          
if size(x,1)==1
    warning(['There is only one sample for ' ...
        'covariance matrix estimation']);
    sig = 0;
    return;
end
            
switch param.method
    case 'mle'
        if isempty(param.weights)
            sig = cov(x,1);
        else
            sig = tz_cov(x,1,param.weights);
        end


        [R,p] = chol(sig);
        if p>0 %if the covariance matrix is not positive definite
            warning(['The covariance matrix is not ' ...
                'positive definite. Calibration is done.']);
            sig = sig+diag(diag(sig)*param.calscale);
        end
        
    case 'ind'
        if isempty(param.weights)
            sig = diag(var(x,1,1));
        else
            for i=1:size(x,2)
                sig(i,i) = ml_wmoment(x(:,i),param.weights);
            end
        end
    case 'idv'
        n = size(x,1);
        k = size(x,2);
        if isempty(param.weights)
            v = sum(var(x,1,1))/k;
        else
            v = ml_wmoment(x(:),repmat(param.weights,k,1));
        end
        sig = diag(zeros(1,k)+v);
    otherwise
        error('Unrecognized covariance estimation method.');
end
