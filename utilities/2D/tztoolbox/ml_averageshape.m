function shape = ml_averageshape(shapes,niter,isshow)
%ML_AVERAGESHAPE Active shape model.
%   SHAPE = ML_AVERAGESHAPE(SHAPES,NITER) returns a structure of shape
%   statistical model learned from the cell array of shapes, SHAPES. Each 
%   element of SHAPES is an array of points. All elements have the same 
%   size. NITER is the number of iteration.
%   
%   SHAPE = ML_AVERAGESHAPE(SHAPES,NITER,0) disables displaying
%   intermediate results.
%
%   See also

%   25-Apr-2005 Initial write T. Zhao
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

if ~exist('isshow','var')
    isshow = 1;
end

baseshape=shapes{1};
tpts{1}=baseshape;
for i=[2:length(shapes)]
    [a,t]=ml_alignshape(shapes{i},baseshape);
    A=[a(1) a(2) 0;-a(2) a(1) 0;[t 1]];
    tpts{i}=tz_deformpts_2d(shapes{i},A);
end

sumshape=zeros(size(shapes{1}));
for i=1:length(shapes)
    sumshape=sumshape+tpts{i};
end

avgshape=sumshape/length(shapes);

for j=1:niter-1
    baseshape=avgshape;

    for i=1:length(shapes)
        [a,t]=ml_alignshape(shapes{i},baseshape);
        A=[a(1) a(2) 0;-a(2) a(1) 0;[t 1]];
        tpts{i}=ml_deformpts_2d(shapes{i},A);
    end
    
    sumshape=zeros(size(shapes{1}));
    for i=1:length(shapes)
        sumshape=sumshape+tpts{i};
    end
    
    avgshape=sumshape/length(shapes);
%     if isshow
%         ml_showpts_2d(avgshape,'ln',1);
%         drawnow
%     end
end

for i=1:length(shapes)
    shapemat(i,:)=tpts{i}(:)'-avgshape(:)';
end

[pccoeff,pcvec]=pca(shapemat);

pcacom=shapemat*pcvec;

shape=struct('avgshape',avgshape,'pcacom',pcacom,'pccoeff',pccoeff, ...
    'pcvec',pcvec,'name','active')
