function [p,refimg,ref,fitdata] = ml_fitdistshape3d(img,param)
%ML_FITDISTSHAPE3D Fit the height of a 3D shape by distances.
%   P = ML_FITDISTSHAPE3D(IMG) returns the parameters of the fitting results.
%   P is a structure with 
%   
%   P = ML_FITDISTSHAPE3D(IMG,PARAM)
%   
%   [P,REFIMG,REF] = ML_FITDISTSHAPE3D(...)
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


if nargin < 1
    error('1 or 2 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,struct('isshow',0,'t1',1.5,'t2',3));

objsize = [];
for i=1:size(img,3)
    objimg = imfill(img(:,:,i),'hole');
    objsize(i) = sum(objimg(:));
    img(:,:,i) = bwperim(objimg);
end

[maxsize,ref] = max(objsize);
refimg = img(:,:,ref);
[hs2,maxdist]=ml_hightcode_3d(img,struct('ref',ref));

hs = hs2;
hs(hs(:,2)<0,:) = [];

fun = struct('funname','ml_projball2','param',param.t1);

x = hs(:,1)/maxdist;
y = hs(:,2);
p.b1 = regress(y,ml_evalfun(x,fun));
p.x = x; p.y = y;
p.t1 = param.t1;
p.fun1 = fun; 
p.fun1.transform.output = ...
    struct('funname','ml_linfun','param',struct('scale',p.b1));

if param.isshow==1
    plot(x,y,'.')
    hold on
    x=0:0.01:1;
    y = p.b1*ml_evalfun(x,fun);
    plot(x,y)
    hold off
end

hs = hs2;
hs(hs(:,2)>0,:) = [];

fun = struct('funname','ml_projball2','param',param.t2);

x = hs(:,1)/maxdist;
y = abs(hs(:,2));
p.b2 = regress(y,ml_evalfun(x,fun));
p.x2 = x; p.y2 = y;
p.t2 = param.t2;
p.fun2 = fun; 
p.fun2.transform.output = ...
    struct('funname','ml_linfun','param',struct('scale',p.b2));

if param.isshow==1
    figure
    plot(x,y,'.')
    hold on
    x=0:0.01:1;
    y = p.b2*ml_evalfun(x,fun);
    plot(x,y)
    hold off
end
