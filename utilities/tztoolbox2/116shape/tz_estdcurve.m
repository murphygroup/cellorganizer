function ds = tz_estdcurve(x,order)
%TZ_ESTDCURVE Estimate derivatives of a curve.
%   DS = TZ_ESTDCURVE(X) returns the derivative of each point in X, which
%   is a 1D array of function values.
%   
%   DS = TZ_ESTDCURVE(X,N) returns N-order derivative of each point in X.
%   
%   See also

%   25-Mar-2005 Initial write TINGZ T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University


if nargin<2
    order=1;
end

if order==1
    ds=[x(2:end),x(1)]-x;
else
%     for i=1:order
%         x=tz_estdcurve(x);
%     end
    ds=tz_estdcurve(tz_estdcurve(x,order-1),1);
end

% switch order
% case 1
%     
% case 2
%     ds=([x(2:end),x(1)]-[x(end),x(1:end-1)])/2;
% end