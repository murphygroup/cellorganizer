function B = tz_multrow(A,x)
%TZ_MULTROW Multiply the a row vector to each row of a matrix.
%   TZ_MULTROW(A,X) returns a matrix in which each row is the array product
%   between teh corresponding row in matrix A and the row vector X.
%   
%   See also ML_ADDROW

%   13-Mar-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

error(tz_genmsg('of','tz_multrow','ml_multrow'));

if nargin < 2
    error('Exactly 2 arguments are required')
end

if(size(A,2)~=size(x,2))
    error('The matrix and the vector must have the same number of columns');
end

B=A.*(ones(size(A,1),1)*x);