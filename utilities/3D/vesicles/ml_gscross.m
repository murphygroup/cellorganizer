function x = ml_gscross(mu1,sigma1,mu2,sigma2)
%ML_GSCROSS The cross points of two gaussian distributions.
%   X = ML_GSCROSS(MU1,SIGMA1,MU2,SIGMA2) returns the values at which two
%   Gaussian distributions have equal density. The two Gaussian
%   distributions are:
%       X1 ~ N(MU1,SIGMA1)
%       X2 ~ N(MU2,SIGMA2)
%
%   See also

%   27-Jan-2006 Initial write T. Zhao
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

if nargin < 4
    error('Exactly 4 arguments are required')
end

c(1) = (1/sigma2^2-1/sigma1^2)/2;
c(2) = -(mu2/sigma2^2-mu1/sigma1^2);
c(3) = (mu2^2/sigma2^2-mu1^2/sigma1^2)/2+log(sigma2/sigma1);

x = roots(c);
