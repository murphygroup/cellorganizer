function y = tz_resetdiag(x,v)
%TZ_RESETDIAG Reset values in the diagonal of a matrix.
%   Y = TZ_RESETDIAG(X) resets the digonal of X to all 0s.
%   
%   Y = TZ_RESETDIAG(X,V) resets the digonal of X to value V. V could be a
%   scalar or vector. If it is a vector, it must have the same size as the
%   diagonal of X.
%   
%   See also

%   23-Aug-2006 Initial write T. Zhao
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


if nargin < 1
    error('1 or 2 arguments are required');
end

if ~exist('v','var')
    v=0;
end

if ~isvector(v)
    error('Wrong data type. The input must be a vector.');
end

if tz_vectype(v)==1
    v = v';
end

y = x-diag(diag(x)-v);


