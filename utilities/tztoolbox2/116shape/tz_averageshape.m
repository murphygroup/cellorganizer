function shape = tz_averageshape(shapes,niter,isshow)
%TZ_AVERAGESHAPE Obsolete. See ML_AVERAGESHAPE.
%   SHAPE = TZ_AVERAGESHAPE(SHAPES,NITER) returns a structure of shape
%   statistical model learned from the cell array of shapes, SHAPES. Each 
%   element of SHAPES is an array of points. All elements have the same 
%   size. NITER is the number of iteration.
%   
%   SHAPE = TZ_AVERAGESHAPE(SHAPES,NITER,0) disables displaying
%   intermediate results.
%
%   See also

%   25-Apr-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

error(tz_genmsg('of','tz_averageshape','ml_averageshape'));

if ~exist('isshow','var')
    isshow = 1;
end

baseshape=shapes{1};
tpts{1}=baseshape;
for i=[2:length(shapes)]
    [a,t]=tz_alignshape(shapes{i},baseshape);
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
        [a,t]=tz_alignshape(shapes{i},baseshape);
        A=[a(1) a(2) 0;-a(2) a(1) 0;[t 1]];
        tpts{i}=tz_deformpts_2d(shapes{i},A);
    end
    
    sumshape=zeros(size(shapes{1}));
    for i=1:length(shapes)
        sumshape=sumshape+tpts{i};
    end
    
    avgshape=sumshape/length(shapes);
    if isshow
        tz_showpts_2d(avgshape,'ln',1);
        drawnow
    end
end

for i=1:length(shapes)
    shapemat(i,:)=tpts{i}(:)'-avgshape(:)';
end

[pccoeff,pcvec]=pca(shapemat);

pcacom=shapemat*pcvec;

shape=struct('avgshape',avgshape,'pcacom',pcacom,'pccoeff',pccoeff, ...
    'pcvec',pcvec,'name','active')
