function xy = tz_cov(x,varargin)
%TZ_COV Obsolete. See ML_COV.
%   COV(X), if X is a vector, returns the variance.  For matrices,
%   where each row is an observation, and each column a variable,
%   COV(X) is the covariance matrix.  DIAG(COV(X)) is a vector of
%   variances for each column, and SQRT(DIAG(COV(X))) is a vector
%   of standard deviations. COV(X,Y), where X and Y are
%   vectors of equal length, is equivalent to COV([X(:) Y(:)]). 
%   
%   COV(X) or COV(X,Y) normalizes by (N-1) where N is the number of
%   observations.  This makes COV(X) the best unbiased estimate of the
%   covariance matrix if the observations are from a normal distribution.
%
%   COV(X,1) or COV(X,Y,1) normalizes by N and produces the second
%   moment matrix of the observations about their mean.  COV(X,Y,0) is
%   the same as COV(X,Y) and COV(X,0) is the same as COV(X).
%
%   The mean is removed from each column before calculating the
%   result.
%
%   See also CORRCOEF, STD, MEAN.

%   J. Little 5-5-86
%   Revised 6-9-88 LS 3-10-94 BJ
%   Copyright 1984-2001 The MathWorks, Inc. 
%   $Revision: 1.3 $  $Date: 2006/09/11 17:30:34 $

error(tz_genmsg('of','tz_cov','ml_cov'));

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

