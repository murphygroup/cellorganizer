function Y=ml_pdist(X,s,t)

%ML_PDIST Pairwise distance between observations.
%Modified from pdist function
%   Y =ML_ PDIST(X,METRIC) returns a vector which contains all the
%   distances between each pair of observations in X computed using
%   the given METRIC.  X is a M by N matrix, treated as M observations
%   of N variables. Since there are M*(M-1)/2 pairs of observations in
%   X, the size of Y is M*(M-1)/2 by 1.  The default metric is
%   'EUCLID'.  The available metrics are:
%
%      'euclid'    --- Euclidean metric
%      'seuclid'   --- Standardized Euclid metric
%      'cityblock' --- City Block metric
%      'mahal'     --- Mahalanobis metric
%      'minkowski' --- Minkowski metric
%
%   Probability distance is a distance defined by tingz for hypothesis test 
%   results. It's not a metric distance. 
%
%   Y = PDIST(X, 'minkowski', p) specifies the exponents in the
%   'Minkowski' computation. When p is not given, p is set to 2.
%
%   The output Y is arranged in the order of ((1,2),(1,3),..., (1,M),
%   (2,3),...(2,M),.....(M-1,M)).  i.e. the upper right triangle of
%   the M by M square matrix. To get the distance between observation
%   i and observation j, either use the formula Y((i-1)*(M-i/2)+j-i)
%   or use the helper function Z = SQUAREFORM(Y), which will return a
%   M by M symmetric square matrix, with Z(i,j) equaling the distance
%   between observation i and observation j.
%
%   See also SQUAREFORM, LINKAGE


if nargin >= 2
   if length(s) > 2
      error('Unrecognized metric');
  else 
      s = lower(s(1:2));
   end
else
   s = 'eu';
end


if nargin < 3
    if all(s=='mi')
        t = 2; 
    end
    
    if all(s=='pv')
        t=0.5;
    end
elseif t <= 0
    error('The third argument has to be positive.');
end


[m, n] = size(X);

if m < 2
   error('The first argument has to be a numerical matrix with at least two rows');
end

p = (m-1):-1:2;
I = zeros(m*(m-1)/2,1);
I(cumsum([1 p])) = 1;
I = cumsum(I);
J = ones(m*(m-1)/2,1);
J(cumsum(p)+1) = 2-p;
J(1)=2;
J = cumsum(J);

Y = (X(I,:)-X(J,:))';
p = [];  % no need for I J p any more.

switch s
case 'eu' % Euclidean
   Y = sum(Y.^2,1);
   Y = sqrt(Y);
case 'se' % Standadized Euclidean
   D = diag(1./var(X));
   Y = sum(D*(Y.^2),1);
   Y = sqrt(Y);
case 'ci' % City Block
   Y = sum(abs(Y),1);
case 'ma' % Mahalanobis
   v = inv(cov(X));
   Y = sqrt(sum((v*Y).*Y,1));   
case 'm2' %Two sample Mahalanobis
    N=size(X,1);
    n1=t;
    n2=N-t;
    cov1=cov(X(1:t,:));
    cov2=cov(X((t+1):end,:));
    s = ((n1-1)*cov1 + (n2-1)*cov2)/(n1+n2-2);
    v = inv(s);
    Y = sqrt(sum((v*Y).*Y,1)); 
case 'mi' % Minkowski
   Y = sum(abs(Y).^t,1).^(1/t);
otherwise
   error('no such method.');
end

I = []; J = []; 
