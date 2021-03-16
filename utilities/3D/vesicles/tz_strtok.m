function ts = tz_strtok(s,d)
%TZ_STRTOK Find tokens in string.
%   TS = TZ_STRTOK(S) returns all tokens in the string S delimited by
%   "white space". TS is a cell array of these tokens.
%   
%   TS = TZ_STRTOK(S,D) returns all tokens delimited by one of the 
%   characters in D. TS is a cell array of these tokens.
%   
%   Example:
%       tz_strtok('Enjoy using this toolbox!') returns 
%       {'Enjoy','using','this','toolbox!'}
%
%       tz_strtok('Enjoy&use this,toolbox!','&,!') returns 
%       {'Enjoy','use this','toolbox'}
%
%   See also STRTOK

%   11-Aug-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University


% Copyright (C) 2007  Murphy Lab
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

if nargin<1
    error('1 or 2 arguments are required');
end

if ~ischar(s)
    error('The fist argument must be a string');
end

if nargin<2
    d=' ';
end

if isempty(d)
    ts{1}=s;
    return
end

ts={};

[token,remain] = strtok(s,d);

while ~isempty(token)
    ts{end+1}=token;
    [token,remain] = strtok(remain,d);
end
    
