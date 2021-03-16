function [mu, coeff, score, latent] = eff_PCA(A)
% Let A be a N*n matrix, N rows represent the observation
% n columns represent the variance.
% Here we calculate the covariance matrix of input A, in order 
% calculate the eigen values and eigen vectors of the A*A'
%   MU - Recentering of the orginal data.
%   COEFF - n*N matrix, correspond to the N normalized eigen vectors of the
%       covariance matrix, i.e., principle components.
%   SCORE - N*N matrix, the projected length on each principle mode. There
%       are N principle modes as maximum.
%   LATENT - N diangonal vector correspond to the eigen values of each
%       principle mode
% 
% The ultimate goal of doing so, is to use the optimized method 
% to compute PCA from A*A' (a N by N matrix) rather than 
% frome A'*A (a n by n matrix). Because the dimension n of the 
% variance is often much bigger than the observation N.
% 
% We can do so by:
% 1. Compute the eigVect and eigVal of the T = (D * D' )
%     Then we have T * eigVect = eigVal * eigVect;
% 2. We can compute the eigVect and eigVal of S = D' * D 
%     by D' * eigVect


% Author: Wei Wang
% September ??, 2010 T. Peng
%
% Copyright (C) 2012 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
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

% Step1 compute the eigVect and eigVal of the T = (D * D' )
[N, n] = size(A);
mu = mean(A);
D = (A - repmat(mu, N,1))';

T = D'*D/N;
[V, eigD] = eig(T);

eigenVect = fliplr(V);
latent = flipud(diag(eigD));

% Step 2 compute A' * eigenVect, in order to calculate the eigenVect
% of the original covariance matrix

coeff = D * eigenVect;
for i = 1:N
    coeff(:,i) = coeff(:,i) / sqrt(latent(i) * N);
end

score = coeff' * D;
