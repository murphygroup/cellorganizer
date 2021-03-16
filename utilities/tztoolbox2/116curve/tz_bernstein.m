function y = tz_bernstein(x,betas)
%TZ_BERNSTEIN Evaluate a polynomial in Bernstein form.
%   Y = TZ_BERNSTEIN(X,BETAS) returns the values of the polynomial at X,
%   which must be a vector. The polynomial has Bernstein form with
%   coeffiecients BETAS, which is a row vector of length N+1, where N is 
%   the order of the polynomial. The function is implemented by the de 
%   Casteljau's algorithm.
%   
%   See also

%   27-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('Exactly 2 arguments are required')
end

nsample = length(x);
if size(x,2)>1
    x = x';
end

n = length(betas)-1;
rbetas{1} = repmat(betas,nsample,1);

for j=1:n
    for i=0:n-j
        rbetas{j+1}(:,i+1) = rbetas{j}(:,i+1).*(1-x)+rbetas{j}(:,i+2).*x;
    end
end

y = rbetas{n+1}(:,1);