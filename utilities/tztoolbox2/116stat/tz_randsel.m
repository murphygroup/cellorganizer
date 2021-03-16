function randsel=tz_randsel(n,m)

%randsel=tz_randsel(n,m)
%OVERVIEW:
%   Select m elements from n integers randomly
%PARAMETERS:
%   n - 1:n
%   m - sel number
%RETURN:
%   randsel -  a vector, sorted
%DESCRIPTION:
%
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   05-NOV-2004 Modified TINGZ
%       - add comments


if n<m
    warning('n should be greater than m')
    return
end

ro=tz_randorder(n);
randsel=sort(ro([1:m]));
