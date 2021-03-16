function s=tz_num2str(x,brk)
%TZ_NUM2STR Convert a numerical vector to a string.
%   TZ_CELL2STR(X) returns a string, which is the comma separated 
%   combination of numbers in the vector X.
%
%   TZ_CELL2STR(X,D) lets a user specify a customized sparator by D.
%
%   Example:
%       tz_cell2str([1 2 3]) returns '1,2,3'   
%       tz_cell2str([1 2 3],'-') returns '1-2-3'
%
%   See also TZ_CELL2STR

%   23-MAY-2003 Initial write TINGZ
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('1 or 2 arguments are required');
end

if all(size(x) > 1)
    error('The 1st argument must be a vector');
end

if isempty(x)
    s=[];
    return
end

if ~exist('brk','var')
    s = tz_cell2str(num2cell(x));
else
    s = tz_cell2str(num2cell(x),brk);
end

% if ~exist('brk','var')
%     brk=',';
% end
% 
% n=length(x);
% s=[];
% for i=1:n-1
%     s=[s,num2str(x(i)),brk];
% end
% s=[s,num2str(x(end))];
