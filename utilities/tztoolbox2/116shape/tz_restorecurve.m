function y = tz_restorecurve(x,order)
%TZ_RESTORECURVE Restore a derived curve.
%   Y = TZ_RESTORECURVE(X,N) returns a curve that has N-order derivative X.
%   
%   See also

%   25-Mar-2005 Initial write  T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end 

y=cumsum(x);

for i=2:order
    y=cumsum(-mean(y)+y);
end
y=[y(end-order+1:end),y(1:end-order)];