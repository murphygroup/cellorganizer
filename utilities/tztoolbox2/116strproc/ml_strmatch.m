function pos = ml_strmatch(strs1,strs2)
%ML_STRMATCH Find a set of strings in anther set of strings.
%   POS = ML_STRMATCH(STRS1,STRS2) returns the indices of where the strings in
%   the cell array STRS1 first occur in the cell array STRS2. If there is no 
%   occurence of a string, its index is 0.
%   
%   See also

%   14-Feb-2007 Initial write T. Zhao
%   Copyright (c) 2007 Murphy Lab
%   Carnegie Mellon University
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation; either version 2 of the License,
%   or (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301, USA.
%   
%   For additional information visit http://murphylab.web.cmu.edu or
%   send email to murphy@cmu.edu


if nargin < 2
    error('Exactly 2 arguments are required');
end

nstr2 = length(strs2);
[tmp,pos] = ismember(strs1,fliplr(strs2));
indices = find(pos>0);
pos(indices)= nstr2 - pos(indices)+1;
