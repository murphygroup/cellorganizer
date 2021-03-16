function [img,obj]=ml_findmainobj_bw(img,n)
%ML_FINDMAINOBJ_BW Find the biggest object in an image.
%   OBJ = ML_FINDMAINOBJ_BW(IMG) returns a 2-column matrix representing 
%   the extracted object from a binary image IMG. 
%   OBJ = TZ_FINDMAINOBJ_BW(IMG,N) specifies the parameter for BWLABEL. 
%   N defaults to 8.

%   30-SEP-2004 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

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
    error('1 or 2 arguments are required')
end

if nargin<2
    n=18;
end

limg=bwlabeln(img,n);
objnum=max(limg(:));
        
lhist=[];
for i=1:objnum
    lhist(i)=sum(sum(sum(limg==i)));
end
[y,maxl]=max(lhist);
img(limg~=maxl)=0;

idx = find(img==1);
[X,Y,Z] = ind2sub(size(img),idx);
obj = [X,Y,Z];
