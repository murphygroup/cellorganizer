function ml_imvis(img,param)
%ML_IMVIS Visualize a 3d image.
%   ML_IMVIS(IMG) visualize the 3d image IMG.
%   
%   ML_IMVIS(IMG,PARAM) allows users to customize the visualization by the
%   structure PARAM, which has the following fields:
%       'isovalue' - value for the isosurface to display. Default: 1.
%       'isedge' - shows the outer boundary of IMG if it is 1. It expects
%           IMG to be binary and 'isovalue' has no use in this case. The
%           default value is 0, which means this option is off.
%       'zoffset' - offset along z axis. Default: 0.
%       'FaceColor' - surface color. Default: 'green'.
%       'EdgeColor' - edge color. Default: 'none'.
%       'transparency' - transparency. Default: 1 (not transparent). See ALPHA
%           for more details.
%       'reducevolume' - reduce factor. See REDUCEVOLUME. Default: empty
%           (no reduce).
%       'shrinkfaces' - shrink faces. See SHRINKFACES. Default: empty (no
%          shrink).
%       
%   See also

%   17-Oct-2006 Initial write T. Zhao
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

param = ml_initparam(param,struct('isovalue',1,'isedge',0, ...
    'zoffset',0, 'reducevolume',[],'shrinkfaces',[],'reducepatch',[],...
    'FaceColor','green','EdgeColor','none','transparency',1));

if param.isedge==1
    param.isovalue = 0;
    img = double(img);
    for i=1:size(img,3)
        objimg = imfill(img(:,:,i),'hole');
        distimg = double(bwdist(bwperim(objimg)));
        distimg(objimg==0) = -distimg(objimg==0);
        img(:,:,i) = distimg;
    end
end

[x,y,z] = meshgrid(1:size(img,2),1:size(img,1),1:size(img,3));
z = z+param.zoffset;
v = img;

if ~isempty(param.reducevolume)
    [x,y,z,v] = reducevolume(x,y,z,v,param.reducevolume);
end

fv = isosurface(x,y,z,v,param.isovalue);

if ~isempty(param.shrinkfaces)
    fv = shrinkfaces(fv,param.shrinkfaces);
end

p = patch(fv);

isonormals(x,y,z,v,p);
set(p,'FaceColor',param.FaceColor','EdgeColor',param.EdgeColor);

if ~isempty(param.reducepatch)
    figure;
    h = axes;
    p = copyobj(p,h);
    p = reducepatch(p,param.reducepatch);
end



daspect([1 1 1])
view(3); axis tight
camlight
lightangle(-40,-30);
lighting gouraud
alpha(param.transparency);

