function y = ml_sumfun(x,param)
%ML_SUMFUN The sum of several functions.
%   Y = ML_SUMFUN(X,PARAM) sum the functions containd in the structure PARAM
%   up. The field 'fs' in PAMAM indidates the set of fuctions. It is a cell
%   array of [general function]s, which should have the same types of input
%   and output. PARAM could also have a field called 'ws', which is a vector
%   of weights. So the combined function will be ws(1)*fs{1}+ws(2)*fs{2}+...
%   If the field does not exist or it is empty no weights will be added. 
%   Otherwise its length should equal the number of functions.
%
%   See also

%   10-Nov-2006 Initial write T. Zhao
%   Copyright (c) 2006 Murphy Lab
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

param = ml_initparam(param,struct('ws',[]));

if isempty(param.ws)
    param.ws = ones(1,length(param.fs));
end

y = ml_evalfun(x,param.fs{1})*param.ws(1);

for i=2:length(param.fs)
    y = y+ml_evalfun(x,param.fs{i})*param.ws(i);
end
    
