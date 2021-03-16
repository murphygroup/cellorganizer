function img2 = ml_bwrecon3d(img,param)
%ML_BWRECON3D Reconstruct a 3D image from a 2D image.
%   IMG2 = ML_BWRECON3D(IMG,PARAM) returns a 3D image that is recontructed
%   from a [binary image] img. PARAM is a structure with the following fields:
%       fun1 - [general function] for top height
%       fun2 - [general function] for bottom height
%       maxheight1 - maximum top height
%       maxheigth2 - minimum top height
%
%   See also

%   18-Oct-2006 Initial write T. Zhao
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


if nargin < 2
    error('Exactly 2 argument is required');
end

param = ml_initparam(param,struct('ratio',1));

refimg = img;
distimg = bwdist(refimg);
maxdist = max(max(distimg(imfill(refimg,'hole')==1)));

%fun = struct('funname','ml_projball2','param',param.t1);

%maxheight = floor(param.b1);
maxheight = param.maxheight1;

% ds = -ml_evalfun((0:maxheight)/param.b1,ml_getinvfun(fun))*maxdist;
ds = -ml_evalfun((0:maxheight),ml_getinvfun(param.fun1))*maxdist;

rimg = [];
for i=1:length(ds)
    rimg(:,:,i) = ml_bwdistshape(refimg,ds(i));
end

%fun = struct('funname','ml_projball2','param',param.t2);

%maxheight = floor(param.b2);
maxheight = param.maxheight2;

% ds = -ml_evalfun((1:maxheight)/param.b2,ml_getinvfun(fun))*maxdist;
ds = -ml_evalfun((1:maxheight),ml_getinvfun(param.fun2))*maxdist;
img2 = [];
for i=1:length(ds)
    img2(:,:,i) = ml_bwdistshape(refimg,ds(i));
end

img2 = flipdim(img2,3);

img2(:,:,end+1:end+size(rimg,3)) = rimg;

