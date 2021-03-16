function y = ml_pcatrans(x,param)
%ML_PCATRANS PCA transformation.
%   Y = ML_PCATRANS(X,PARAM) returns the PCA transformation of X. PARAM is a
%   structure with the following fields:
%       'basevec' - base vectors. Each column is a base vector.
%       'offset' - offset of X before transformation. The offset will be
%            subtracted from X. If it is empty or does 
%            not exist, no offset will be done.
%   
%   See also

%   02-Nov-2006 Initial write T. Zhao
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

param = ml_initparam(param,struct('offset',[]));
if ~isempty(param.offset)
    y = ml_addrow(x,-param.offset);
end

y = y*param.basevec;
