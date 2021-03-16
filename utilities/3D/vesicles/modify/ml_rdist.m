function rdist = ml_rdist(nucEdge,cellEdge,nucCentre,nucAngle,param)
%ML_RDIST returns sampled distance ratio between nuclear edge radius and
% cell edge radius
%   See also

%   10-Jan-2006 Initial write T. Peng
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

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,struct('Nxy',360,'Nz',180));

nucEdge = ml_addrow(nucEdge,-nucCentre);
cellEdge = ml_addrow(cellEdge,-nucCentre);

[nucTh,nucPhi,nucR] = cart2sph(nucEdge(:,1),nucEdge(:,2),nucEdge(:,3)*5);
[cellTh,cellPhi,cellR] = cart2sph(cellEdge(:,1),cellEdge(:,2),cellEdge(:,3)*5);
clear nucEdge
clear cellEdge
nucTh = nucTh - nucAngle;
nucTh(nucTh<=-pi) = nucTh(nucTh<=-pi) + 2*pi;
cellTh = cellTh - nucAngle;
cellTh(cellTh<=-pi) = cellTh(cellTh<=-pi) + 2*pi;

% Sampling
Nz = param.Nz;
Nxy = param.Nxy;
Phi = (1 + cos(-pi/2:pi/Nz:pi/2)) .* (-pi/2:pi/Nz:pi/2);
Phi(1) = []; Phi(end) = [];
Th = (-(Nxy/2-1):(Nxy/2)) * (2*pi/Nxy);
[TH,PHI] = meshgrid(Th,Phi);
% TH = TH(:); PHI = PHI(:);
NR = griddata(nucTh,nucPhi,nucR,TH,PHI);
CR = griddata(cellTh,cellPhi,cellR,TH,PHI);
rdist = CR ./ NR;


