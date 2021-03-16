function y = ml_shift(x,offset)
%ML_SHIFT Circular shift of a matrix.
%   Y = ML_SHIFT(X,OFFSET) return a circularly shifted matrix of X. OFFSET is a
%   scalar of 1x2 vector. If it is a scalar, both dimensions will be shifted by
%   OFFSET. If it is a 1x2 vector, OFFSET(1) is for rows and OFFSET(2) is for 
%   columns. The direction is left and down.
%   
%   See also

%   10-Sep-2006 Initial write T. Zhao
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

if length(offset)==1
    offset = [offset offset];
end

y = circshift(x,offset);
