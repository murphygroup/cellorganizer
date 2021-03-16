function answer = is_ome_tiff( filename )

% Copyright (C) 2018  Murphy Lab
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

answer = false;

try
	tiff=Tiff(filename);
	ImageDescription=tiff.getTag('ImageDescription');
catch
	warning('Fail to get ImageDescription.');
	return
end

if isempty(strfind(ImageDescription,'</OME>'))
	warning('The input tiff image does not has <OME> tag.');
	return
end

%%%
data = bfopen(filename);
omeMeta = data{1,4}; %gets the metadata block
try
	dimX = omeMeta.getPixelsSizeX(0).getValue();
catch
	warning('X dimension does not exist.');
	return
end

try
	dimY = omeMeta.getPixelsSizeY(0).getValue();
catch
	warning('Y dimension does not exist.');
	return
end

try
	dimZ = omeMeta.getPixelsSizeZ(0).getValue();
catch
	warning('Z dimension does not exist.');
	return
end

try
	dimC = omeMeta.getPixelsSizeC(0).getValue();
catch
	warning('C dimension does not exist.');
	return
end

try
	dimT = omeMeta.getPixelsSizeT(0).getValue();
catch
	warning('T dimension does not exist.');
	return
end

if dimX < 0
	warning('x is negative.');
	return
end
if dimY < 0
	warning('y is negative.');
	return
end
if dimZ < 0
	warning('z is negative.');
	return
end
if dimC < 0
	warning('c is negative.');
	return
end
if dimT < 0
	warning('t is negative.');
	return
end

%%%
filename = char(filename);
[~, name, ext1] = fileparts( filename );
[~, ~, ext2] = fileparts( name);
if strcmp(ext2, '.ome')
    % standard TIFF format: .ome.tif, .ome.tiff
    % BigTiff specific extension: .ome.tf2, .ome.tf8, .ome.btf
    if strcmp(ext1, '.tif') || strcmp(ext1, '.tiff') || ...
	    strcmp(ext1, '.tf2') || strcmp(ext1, '.tf8') || ...
	    strcmp(ext1, '.btf')
	answer = true;
    end
end
end

