function map = tp_cartimg(surfmap)
% TP_SURFMAP convert a 3D surface into a 2D map using Cartesian-Cylindrical
% coordinate transformation. DELTA is the resampling step size of theta

% Author: Tao Peng
% Edited: Ivan E. Cao-Berg
%
% Copyright (C) 2011 Murphy Lab
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

[H,P] = size(surfmap);
delta = 2*pi/(P-1);

Phi = -pi:delta:pi;
N = length(Phi);
R = zeros(S,N);

for k = 2:S+1
    % Resampling
    [boundX,boundY] = pol2cart(Phi,surfmap(k,:));
    boundX = round(boundX+ctrX);
    boundY = round(boundY+ctrY);
    boundX(isnan(boundX)) = [];
    boundY(isnan(boundY)) = [];
    for i = 1:length(boundX)
        recons([boundX(i),boundY(i)],k) = 1;
    end
end


[H,W,S] = size(img);
