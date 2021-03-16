function v=tz_mat2vec(A,bycol)
%TZ_MAT2VEC Convert a matrix to a vector.
%   TZ_MAT2VEC(A) is the same as A(:), which reshape the matrix A to a
%   vector columnwise.
%   
%   TZ_MAT2VEC(A,BYCOL) uses BYCOL to specify columnwise or rowwise
%   convertion. If BYCOL is 0, then the vector will have rowwise order.
%   Otherwise, it is the same as TZ_MAT2VEC(A).

%   ??-???-???? Initial write T. Zhao
%   07-NOV-2004 Modified T. Zhao
%       - add parameters bycol
%   Copyright (c) Murphy Lab, Carnegie Mellon University


if nargin < 1
    error('1 or 2 arguments are required')
end

if ~exist('bycol','var')
    bycol=1;
end

if bycol==0
    A=A';
end

v=A(:);