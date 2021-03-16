function xy = ml_cov(x,varargin)
%ML_COV Covariance matrix. (Modified from Matlab COV funcion)
%   ML_COV(X), if X is a vector, returns the variance.  For matrices,
%   where each row is an observation, and each column a variable,
%   ML_COV(X) is the covariance matrix.  DIAG(ML_COV(X)) is a vector of
%   variances for each column, and SQRT(DIAG(ML_COV(X))) is a vector
%   of standard deviations.
%   
%   ML_COV(X) onormalizes by (N-1) where N is the number of
%   observations.  This makes ML_COV(X) the best unbiased estimate of the
%   covariance matrix if the observations are from a normal distribution.
%
%   ML_COV(X,1) normalizes by N and produces the second
%   moment matrix of the observations about their mean.  ML_COV(X,Y,0) is
%   the same as ML_COV(X,Y) and ML_COV(X,0) is the same as ML_COV(X).
%
%   ML_COV(X,FLAG,WEIGHTS) calculate weighted covariance matrix.
%
%   ML_COV(X,FLAG,WEIGHTS,MU) calculate the covariance matrix based on the
%   predefined mean MU.
%
%   The mean is removed from each column before calculating the
%   result.
%
%   See also CORRCOEF, STD, MEAN.

%   J. Little 5-5-86
%   Revised 6-9-88 LS 3-10-94 BJ
%   Copyright 1984-2001 The MathWorks, Inc. 
%   $Revision: 1.2 $  $Date: 2006/11/17 19:11:38 $

% Copyright (C) 2007  Murphy Lab
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

if nargin==0, error('Not enough input arguments.'); end
if nargin>4, error('Too many input arguments.'); end
if ndims(x)>2, error('Inputs must be 2-D.'); end

nin = nargin;
if nin>2
    flag=varargin{1};
else
    flag=0;
end

if nin==3
    weights=varargin{2};
else
    weights=[];
end


if length(x)==prod(size(x))
    x = x(:);
end

[m,n] = size(x);

if m==1,  % Handle special case
    xy = 0;
else
    if ~isempty(weights)
        wx=repmat(weights,1,size(x,2)).*x*m/sum(weights);
    else
        wx=x;
    end
    if nin==4
        xbar = varargin{3};
    else
        xbar = [];
    end
    if isempty(xbar)
        xbar=mean(wx,1);
    end
    
    xc = wx - repmat(xbar,m,1);  % Remove mean
%     if ~isempty(weights)
%         wxc=diag(weights)*xc;
%     else
    wxc=xc;
%     end
    if flag
        xy = wxc' * xc / m;
    else
        xy = wxc' * xc / (m-1);
    end
end

