function nucimg = ml_mapnuctex(texModel,nucEdge)
% ML_MAPNUCTEX Synthesize nuclear texture and map it into nuclear boundary.
%   IMG = ML_GENTEX(MODEL) returns a nuclear image with texture that is
%   synthesize from the texture model TEXMODEL filled in the nuclear
%   boundary NUCEDGE.
%   
%   See also ML_TEXMODEL
%
%   22-Sep-2008 Initial write T. Peng
%   Copyright (c) 2006 Murphy Lab
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

% Find the size of the nuc
[Row,Col] = find(nucEdge);
top = min(Row);
bottom = max(Row);
left = min(Col);
right = max(Col);
nucSize = [bottom-top+1,right-left+1];
nucSize = ceil(nucSize/64)*64;

% Synthesize texture image
param = struct('size',nucSize);
teximg = ml_gentex(texModel,param);

% Map the texture
nucimg = zeros(size(nucEdge));
nucimg(top:top+nucSize(1)-1,left:left+nucSize(2)-1) = teximg;
nucForeground = imfill(nucEdge,'hole');
nucimg(~nucForeground) = 0;