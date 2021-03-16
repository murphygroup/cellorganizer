function y = tz_sigmoid(x,threshold,slope)
%TZ_SIGMOID Sigmoid function.
%   Y = TZ_SIGMOID(X) returns 1/(1+e^(-X)).
%   
%   Y = TZ_SIGMOID(X,THRESHOLD) returns 1/(1+e^(threshold-X)).
%
%   Y = TZ_SIGMOID(X,THRESHOLD,SLOPE) returns 1/(1+e^((threshold-X)*SLOPE)).

%   09-Oct-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('1-3 arguments are required')
end

if nargin < 2
    threshold = 0;
end

if nargin < 3
    slope = 1;
end

y = 1./(1+exp(slope*(threshold-x)));