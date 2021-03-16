function [knots,coefs] = ml_sp2feat(sp)
%ML_SP2FEAT Convert B-spline structure into features.
%   KNOTS = ML_SP2FEAT(SP) returns the knots of the spline SP.
%   
%   [KNOTS,COEFS] = ML_SP2FEAT(...) also returns coeffiecients. Both KNOTS
%   and COEFS are scalars or row vectors. KNOTS could also be empty.
%   
%   See also

%   28-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

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

if nargin < 1
    error('Exactly 1 argument is required')
end

if sp.number-sp.order>0
    knots = sp.knots(sp.order+1:sp.number);
else
    knots = [];
end

coefs = sp.coefs;
