function seg = activcon_seg(image,startingmask,display, alpha, maxiter)
%ACTIVCON_SEG Region Based Active Contours

% Author: Shawn Lankton
% 
% Copyright (C) 2012 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% March 15, 2012 R.F. Murphy Use bounding box to speed up calcs
% March 21, 2012 R.F. Murphy Add display argument
% July 1, 2012 G.R. Johnson Allows 3D images
% April 23, 2013 D. Sullivan Added iteration option to the 3D segmentation,
% default 1000
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

if ~exist('alpha', 'var') || isempty(alpha)
    alpha = 0.7;
end

%D. Sullivan 4/23/13 added maxiter variable with default of 1000
if ~exist('maxiter', 'var') || isempty(maxiter)
    maxiter = 1000;
end

I =double(image);
seg = zeros(size(I));
[X,Y]=find(I>0);
border = 10;
% find bounding box but leave some black space around it if possible
Xbox = [max(1,min(X)-border):min(size(I,1),max(X)+border)];
Ybox = [max(1,min(Y)-border):min(size(I,2),max(Y)+border)];

if length(Xbox)<2*border | length(Ybox)<2*border
    return
end

if ~exist('startingmask','var') | isempty(startingmask) | max(max(startingmask))==0
    startingmask = zeros(size(I));          %-- create initial mask

    startingmask(Xbox(border):Xbox(end-border),Ybox(border):Ybox(end-border),:) = 1;
end

%-- Run segmentation
%D. Sullivan 4/23/13 Changed call to region_seg to allow for user specified
%maxiter rather than always using 1000. default still 1000
% seg(Xbox,Ybox,:) = region_seg(I(Xbox,Ybox,:), startingmask(Xbox,Ybox,:), 1000, alpha, display, 0); 
seg(Xbox,Ybox,:) = region_seg(I(Xbox,Ybox,:), startingmask(Xbox,Ybox,:), maxiter, alpha, display, 0); 
