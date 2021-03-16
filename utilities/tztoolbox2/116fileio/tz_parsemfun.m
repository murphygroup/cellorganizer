function s = tz_parsemfun(filename)
%TZ_PARSEMFUN Parse matlab function.
%   S = TZ_PARSEMFUN(FILENAME) return a structure S, which has three
%   fields:
%       inputs: a cell array of input argument names
%       outputs: a cell array of output argment names
%       funname: function name

%   15-Apr-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argumentis required')
end

fullpath=which(filename);
funlines = textread(fullpath,'%s','delimiter','\n','whitespace','');
s=tz_parsemfunlines(funlines);