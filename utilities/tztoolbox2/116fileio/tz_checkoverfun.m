function s = tz_checkoverfun(funname)
%TZ_CHECKOVERFUN Check overloaded functions.
%   S = TZ_CHECKOVERFUN(FUNNAME) checks the overloaded functions with the
%   name FUNNAME. Only functions that can be found in the current paths 
%   will be checked. The returned value S is a structure and has the 
%   following fields:
%       'nfun' - number of functions with the name FUNNAME
%       'funlist' - a [string array] of all found functions. Each element
%           is the fullpath of the corresponding funciton.
%       'diff' - a binary matrix indicating whether the functions are the
%           same. If ith row and jth col is 1, then S.funlist{i} and
%           S.funlist{j} are different. Otherwise they are the same.
%
%   See also

%   02-Mar-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

s.funlist = which([funname '.m'],'-all');

if isempty(s.funlist)
    warning(['No such a function with the name ' funname]);
    return;
end

s.nfun = length(s.funlist);
s.diff = zeros(s.nfun,s.nfun);

for i=1:s.nfun
    for j=i+1:s.nfun
       [s.diff(i,j),msg] = unix(['diff ' s.funlist{i} ' ' s.funlist{j}]);
    end
end

%Make it symmetric
s.diff = s.diff+s.diff';