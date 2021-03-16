function obj = ml_gaussrndobj(sigma,N)
%ML_GAUSSOBJ An object from Gaussian distribution.
%   OBJ = ML_GAUSSOBJ(SIGMA,N) returns an object that is sampled from a
%   Gaussian distribution which has covariance matrix SIGMA and total
%   number of samples N.
%   
%   See also

%   26-Jan-2010 Initial write T. Peng
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

if nargin < 2
    error('Exactly 2 argument is required')
end

N = round(N);
samples = mvnrnd(zeros(N,size(sigma,1)),sigma);
samples = round(samples);
[coords,I,J] = unique(samples,'rows');
intensity = hist(J,max(J));
obj = [coords,intensity'];
