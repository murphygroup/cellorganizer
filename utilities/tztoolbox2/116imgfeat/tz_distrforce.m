function fs=tz_distrforce(X,option)
%TZ_DISTRFORCE Cacluate connecting forces for an image.
%   FS = TZ_DISTRFORCE(X,OPTION) retures a vector of connecting forces
%   between adjacent rows of X. The last element of FS is the force 
%   between the last row and the first row. See TZ_CONNNECTFORCE for
%   details about OPTION.

%   ??-???-???? Initial write T. Zhao
%   04-NOV-2004 Modified T. Zhao
%       - add comments
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

if size(X,1)<=1
    warning('a line can not be cut')
end

for i=1:size(X,1)-1
    fs(i)=tz_connectforce(X(i,:),X(i+1,:),option);
end

fs(end+1)=tz_connectforce(X(end,:),X(1,:),option);


