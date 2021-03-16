% ML Feature Calculation
%  
% Top Level files:
%   ml_featset - preprocess image and call ml_features
%   ml_features - call the feature calculation code for each feature set
%  
% Morphological Features
%   ml_imgfeatures - calculate the object related features
% Geometric Features
%   ml_hullfeatures - Calculate hull features
%   ml_imgconvhull -  returns the convex hull of BWIMAGE as an image
%   ml_setimgptspixel - sets values for specified pixels
%   ml_getlinepts - get all discrete points along a line
% Edge Features
%   ml_imgedgefeatures - calculate edge features
% Skeleton Features
%   ml_imgskelfeats - calculate skelenton features
%   ml_objskelfeats  -   Calculate skeleton features for the object OBJIMG.
%   ml_mmthin - rewrite the mmthin function in the morphological toolbox
%   ml_find_branch_points - find the branch points of the skeleton
% Zernike Moment Features
%   ml_zernike     - Calculate Zernike Moment Features
%   ml_imgmoments  - calculates the moment MXY for IMAGE
% Wavelet Features
%   ml_wavefeatures - Calculate Wavelet Features
% Garbor Features
%   ml_gaborfeat - calculate Garbor features
% Others:
%   ml_imgcentmoments - calculates the central moment MUxy for IMAGE

% Copyright (C) 2006  Murphy Lab
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


