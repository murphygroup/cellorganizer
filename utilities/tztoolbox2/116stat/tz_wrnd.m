function alphas=tz_wrnd(k,n)
%TZ_WRND Generates random weights for a mixture model.
%   ALPHAS = TZ_WRND(K,N) returns a NxK coefficient matrix in which each 
%   row is a sample of mixture model weights. The number of components is
%   K.
%   
%   See also

%   26-MAY-2004 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

for i=1:n
    s=unifrnd(0,1,1,k);
    alphas(i,:)=s/sum(s);
end