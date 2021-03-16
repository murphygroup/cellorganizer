function [vs,a] = ml_covrot(sigma)
%ML_COVROT Covrance matrix rotation.
%   VS = ML_COVROT(SIGMA) returns the variances along the major and minor axis.
%   
%   [VS,A] = ML_COVROT(...) also returns the angle from which SIGMA is rotated 
%   to VS. It has the unit radian.
%   
%   See also

%   06-Feb-2007 Initial write T. Zhao
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


if nargin < 1
    error('1 or 2 arguments are required');
end

[u,s,v] = svd(sigma);

vs = diag(s);
a = acos(u(1,1));


