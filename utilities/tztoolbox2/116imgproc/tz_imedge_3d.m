function edgeimg = tz_imedge_3d(img,param)
%TZ_IMEDGE_3D
%   EDGEIMG = TZ_IMEDGE_3D(IMG)
%   
%   EDGEIMG = TZ_IMEDGE_3D(IMG,PARAM)
%   
%   See also

%   16-Oct-2006 Initial write T. Zhao
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


if nargin < 1
    error('1 or 2 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

img = ml_3dpreprocess(img,[],'none','yesbgsub','nih');

img = img>0;

for i=1:size(img,3)
    edgeimg(:,:,i) = bwperim(img(:,:,i));
end
