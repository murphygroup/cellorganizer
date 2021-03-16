function B=tz_addcol(A,x)
%TZ_ADDCOL add the same vector to each column of a matrix
%   B=TZ_ADDCOL(A,X) adds the vector X to each row of A and 
%   returns the sum. A and X must have the same number of 
%   columns.

if(size(A,1)~=size(x,1))
    error('The matrix and the vector must have the same number of rows');
end

B = ml_addrow(A',x');
B = B';
