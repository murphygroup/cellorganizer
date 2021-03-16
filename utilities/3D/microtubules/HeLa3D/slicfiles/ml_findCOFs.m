function [totalCOF, objCOFs] = ml_findCOFs( objlist, img)

% [TOTALCOF, OBJCOFS] = ML_FINDCOFS( OBJLIST, IMG)
%
% Finds the Center Of Fluorescence (COF) for each object in the
% OBJLIST list of objects, as well as the overall COF.
% IMG is the graylevel image (prior to binarization) that the 
% objects were found from. This is necessary because we want
% to weight the pixel coordinates with their respective graylevel
% values.

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

Size = size(img);
NObj = length( objlist);
objCOFs = zeros( 3, NObj);
TotalSum = [0;0;0];
TotalGray = 0;
for o = 1 : NObj % For each object
    % Get the (y,x,z) coordinates
    v = double( objlist{o}.voxels);
    % Get the graylevel values
    idx = sub2ind( Size, v(1,:), v(2,:), v(3,:));
    graylevel = double( img( idx));
    GraySum = sum( graylevel);
    % Do graylevel weighted coordinate sum
    graylevel = repmat( graylevel, [3 1]);
    CoordSum = sum(v .* graylevel, 2);
    % Find object COFs
    objCOFs(:,o) = CoordSum/GraySum;
    % Add up for Total COF
    TotalSum = TotalSum + CoordSum;
    TotalGray = TotalGray + GraySum;
end
totalCOF = TotalSum/TotalGray;
