function ds = tz_pdist2(x1,x2)
%TZ_PDIST2 Calcuate pairwise Euclidean distances of two matrices.
%   DS = TZ_PDIST2(X1,X2) returns a matrix, in which the ith and jth row is
%   the Euclidean distance between the ith row of X1 and the jth row of X2.
%   X1 and X2 are [feature matrix](ce)s and they should have the same 
%   See also

%   03-Mar-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('Exactly 2 arguments are required')
end

x1square = sum(x1.^2,2);
x2square = sum(x2.^2,2);
crossprod = x1*x2';

ds = sqrt(repmat(x1square,[1,size(x2,1)])+ ...
    repmat(x2square',[size(x1,1),1])-crossprod*2);

