function randorder = tz_randorder(n)

%function randorder = tz_randorder(n)
%
%OVERVIEW:
%   The permutation of consecutive numbers.
%PARAMETERS:
%   n - a random permutation of the integers from 1 to n
%RETURN:
%   random numbers
%DESCRIPTION:
%   the same as randperm
%
%HISTORY:
%   ??-???-???? Initial write TINGZ

error('Function tz_randorder is out of date. Please use randperm');
[b,randorder]=sort(rand(1,n),2);
