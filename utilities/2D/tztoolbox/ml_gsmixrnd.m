function R=ml_gsmixrnd(ps,mus,covs,mm)
%ML_GSMIXRND Sampling from gaussian mixture.
%   R = ML_GSMIXRND(PS,MUS,COVS,MM) returns data sampled from a gaussian
%   mixture distribution, which is determined by PS, MUS and COVS. PS is
%   a vector the coefficients of the compoents and has the length is the
%   same as the number of components. See TZ_MNORND for more details. MS
%   is the a matrix of means of the gaussian distribution. Each column of
%   MS is the mean for a guassian. COVS is a 3D matrix of covariances for
%   the components. COVS(:,:,K) is for the Kth component. The generated
%   data has MM rows and each row is one sample.
%   
%   See also ML_MNORND

%   ??-???-???? Initial write T. Zhao
%   04-NOV-2004 Modified T. Zhao
%       - add comments
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

if nargin < 4
    error('Exactly 4 arguments are required')
end

nums=ml_mnornd(mm,ps,1);
R=[];
for i=1:length(nums)
    if nums(i)>0
        R=[R;mvnrnd(mus(:,i),covs(:,:,i),nums(i))];
    end
end