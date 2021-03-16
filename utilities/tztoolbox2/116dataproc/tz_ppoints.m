function x = tz_ppoints(x1,x2,n)
%TZ_PPOINTS Generates points with equal interval
%   X = TZ_PPOINTS(X1,X2,N) return a row vector of N points betwee X1 and
%   X2.
%   
%   See also

%   02-Nov-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('Exactly 3 arguments are required')
end

if n==1
    x = x1;
    return;
end

if n==2
    x = [x1 x2];
    return;
end

if x1==x2
    x = zeros(1,n)+x1;
    return;
end

interval = (x2-x1)/(n-1);
x = x1:interval:x2;
x(end) = x2;