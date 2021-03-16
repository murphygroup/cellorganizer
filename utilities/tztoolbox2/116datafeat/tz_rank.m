function index=tz_rank(A)
%TZ_RANK Rank a matrix along rows.
%   INDEX = TZ_RANK(A) returns a ranking matrix of A. Each row of INDEX is
%   the rank of the corresponding row in A.
%   
%   Example:
%       TZ_RANK([1 3 6 8 5;1 3 5 2 4]) returns [1 2 4 5 3;1 3 5 2 4].

%   ??-???-???? Initial write T. Zhao
%   05-NOV-2004 Modfied T. Zhao
%       - add comments
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

[m,n]=size(A);
[B,idx]=sort(A,2);
for i=1:m
    index(i,idx(i,:))=1:n;
end