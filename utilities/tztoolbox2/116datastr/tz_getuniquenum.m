function y = tz_getuniquenum(x)
%TZ_GETUNIQUENUM Find unique elements in an integer vector.
%   Y = TZ_GETUNIQUENUM(X) returns a vector which contains all unique
%   elements in the vector X, which must be a vector of natural numbers.
%   This function can be replace by UNIQUE. But TZ_GETUNIQUENUM is faster.

%   31-May-2005  Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University


if nargin < 1
    error('Exactly 1 argument is required')
end

nonsel=zeros(1,max(x));
nonsel(x)=1;
y=find(nonsel==1);