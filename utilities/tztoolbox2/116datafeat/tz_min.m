function [c,m,n] = tz_min(A)
%TZ_MIN Obsolete. See ML_MIN.
%   C = TZ_MIN(A) returns the minimal value in the matrix A.
%   
%   [C,M,N] = TZ_MIN(...) also returns the row and column indices of the
%   minimal value.

%   4-NOV-2003  Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

error(tz_genmsg('of','tz_min','ml_min'));

[c,I]=min(A);
[c,J]=min(c);
m=I(J);
n=J;
